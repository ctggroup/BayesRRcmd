#include "limitsequencegraph.hpp"

#include "BayesRBase.hpp"
#include "compression.h"
#include "marker.h"
#include "markerbuilder.h"

#include <iostream>

LimitSequenceGraph::LimitSequenceGraph(size_t maxParallel)
    : AnalysisGraph(maxParallel)
    , m_graph(new graph)
{
    // Decompress the column for this marker
    auto f = [this] (Message msg) -> Message {
        if (m_bayes->preloaded()) {
            msg.marker = m_bayes->marker(msg.snp);
        } else {
            std::unique_ptr<MarkerBuilder> builder{m_bayes->markerBuilder()};
            builder->initialise(msg.snp, msg.numInds);
            const auto index = m_bayes->indexEntry(msg.snp);
            if (m_bayes->compressed()) {
                builder->decompress(m_bayes->compressedData(), index);
            } else {
                builder->read(m_bayes->preprocessedFile(), index);
            }
            msg.marker.reset(builder->build());
        }
        return msg;
    };
    // Do the decompression work on up to maxParallel threads at once
    m_decompressNode.reset(new function_node<Message, Message>(*m_graph, m_maxParallel, f));

    // The sequencer node enforces the correct ordering based upon the message id
    m_ordering.reset(new sequencer_node<Message>(*m_graph, [] (const Message& msg) -> unsigned int {
        return msg.id;
    }));

    m_ordering2.reset(new sequencer_node<Message>(*m_graph, [] (const Message& msg) -> unsigned int {
        return msg.id;
    }));

    // Do not allow predecessors to carry on blindly until later parts of
    // the graph have finished and freed up some resources.
    m_limit.reset(new limiter_node<Message>(*m_graph, m_maxParallel));

    auto g = [this] (Message msg) -> continue_msg {
        // Delegate the processing of this column to the algorithm class
        m_bayes->processColumn(msg.marker.get());

        // Signal for next decompression task to continue
        return continue_msg();
    };
    // The sampling node is enforced to behave in a serial manner to ensure that the resulting chain
    // is ergodic.
    m_samplingNode.reset(new function_node<Message>(*m_graph, serial, g));

    // Set up the graph topology:
    //
    // orderingNode -> limitNode -> decompressionNode (parallel) -> orderingNode -> samplingNode (sequential)
    //                      ^                                                           |
    //                      |___________________________________________________________|
    //
    // This ensures that we always run the samplingNode on the correct order of markers and signal back to
    // the parallel decompression to keep it constantly fed. This should be a self-balancing graph.
    make_edge(*m_ordering, *m_limit);
    make_edge(*m_limit, *m_decompressNode);
    make_edge(*m_decompressNode, *m_ordering2);
    make_edge(*m_ordering2, *m_samplingNode);

    // Feedback that we can now decompress another column
    make_edge(*m_samplingNode, m_limit->decrement);
}

void LimitSequenceGraph::exec(BayesRBase *bayes,
                              unsigned int numInds,
                              unsigned int numSnps,
                              const std::vector<unsigned int> &markerIndices)
{
    if (!bayes) {
        std::cerr << "Cannot run LimitSequenceGraph without bayes" << std::endl;
        return;
    }

    // Set our Bayes for this run
    m_bayes = bayes;

    // Reset the graph from the previous iteration. This resets the sequencer node current index etc.
    m_graph->reset();

    // Push some messages into the top of the graph to be processed - representing the column indices
    for (unsigned int i = 0; i < numSnps; ++i) {
        Message msg = { i, markerIndices[i], numInds };
        m_ordering->try_put(msg);
    }

    // Wait for the graph to complete
    m_graph->wait_for_all();

    // Clean up
    m_bayes = nullptr;
}
