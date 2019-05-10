#include "parallelgraph.h"

#include "compression.h"
#include "BayesRBase.hpp"
#include "marker.h"
#include "markerbuilder.h"

#include <iostream>

ParallelGraph::ParallelGraph(size_t maxParallel)
    : AnalysisGraph(maxParallel)
    , m_graph(new graph)
{
    // Decompress the column for this marker then process the column using the algorithm class
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

        // Delegate the processing of this column to the algorithm class
        const auto betas = m_bayes->processColumnAsync(msg.marker.get());

        msg.old_beta = std::get<0>(betas);
        msg.beta = std::get<1>(betas);
        msg.deltaEps= std::get<2>(betas);//we update the deltaEps of the message
        return msg;
    };

    // Do the decompress and compute work on up to maxParallel threads at once
    m_asyncComputeNode.reset(new function_node<Message, Message>(*m_graph, m_maxParallel, f));

    // Decide whether to continue calculations or discard
    auto g = [] (decision_node::input_type input,
                 decision_node::output_ports_type &outputPorts) {

        std::get<0>(outputPorts).try_put(continue_msg());

        if (input.old_beta != 0.0 || input.beta != 0.0) {
            // Do global computation
            std::get<1>(outputPorts).try_put(std::move(input));
        } else {
            // Discard
            std::get<0>(outputPorts).try_put(continue_msg());
        }
    };

    m_decisionNode.reset(new decision_node(*m_graph, m_maxParallel, g));

    // Do global computation
    auto h = [this] (Message msg) -> continue_msg {
      m_bayes->updateGlobal(msg.marker.get(), msg.old_beta, msg.beta,msg.deltaEps);//now we pass deltaEps to the global
        return continue_msg();
    };
    // Use the serial policy
    m_globalComputeNode.reset(new function_node<Message>(*m_graph, serial, h));

    // Limit the number of parallel computations
    m_limit.reset(new limiter_node<Message>(*m_graph, m_maxParallel));

    // Enforce the correct order, based on the message id
    m_ordering.reset(new sequencer_node<Message>(*m_graph, [] (const Message& msg) -> unsigned int {
        return msg.id;
    }));

    // Set up the graph topology:
    //
    // orderingNode -> limitNode -> decompressionAndSamplingNode (parallel)
    //                      ^                   |
    //                      |___discard____decisionNode (parallel)
    //                      ^                   |
    //                      |                   | keep
    //                      |                   |
    //                      |______________globalCompute (serial)
    //
    // Run the decompressionAndSampling node in the correct order, but do not
    // wait for the most up-to-date data.
    make_edge(*m_ordering, *m_limit);
    make_edge(*m_limit, *m_asyncComputeNode);

    // Feedback that we can now decompress another column, OR
    make_edge(*m_asyncComputeNode, *m_decisionNode);
    make_edge(output_port<0>(*m_decisionNode), m_limit->decrement);
    // Do the global computation
    make_edge(output_port<1>(*m_decisionNode), *m_globalComputeNode);

    // Feedback that we can now decompress another column
    make_edge(*m_globalComputeNode, m_limit->decrement);
}

void ParallelGraph::exec(BayesRBase *bayes,
                              unsigned int numInds,
                              unsigned int numSnps,
                              const std::vector<unsigned int> &markerIndices)
{
    if (!bayes) {
        std::cerr << "Cannot run ParallelGraph without bayes" << std::endl;
        return;
    }

    // Set our Bayes for this run
    m_bayes = bayes;

    // Do not allow Eigen to parallalize during ParallelGraph execution.
    const auto eigenThreadCount = Eigen::nbThreads();
    Eigen::setNbThreads(0);

    // Reset the graph from the previous iteration. This resets the sequencer node current index etc.
    m_graph->reset();

    // Push some messages into the top of the graph to be processed - representing the column indices
    for (unsigned int i = 0; i < numSnps; ++i) {
        Message msg = { i, markerIndices[i], numInds };
        m_ordering->try_put(msg);
    }

    // Wait for the graph to complete
    m_graph->wait_for_all();

    // Turn Eigen threading back on.
    Eigen::setNbThreads(eigenThreadCount);

    // Clean up
    m_bayes = nullptr;
}
