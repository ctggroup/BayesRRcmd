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
      // Decompress the column
      std::unique_ptr<MarkerBuilder> builder{m_bayes->markerBuilder()};
      builder->initialise(msg.snp, msg.numInds);
      const auto index = m_bayes->indexEntry(msg.snp);
      if (m_bayes->compressed()) {
            builder->decompress(m_bayes->compressedData(), index);
        } else {
            builder->read(m_bayes->preprocessedFile(), index);
          }
      msg.marker.reset(builder->build());
      return msg;
      };

    // Do the decompression work on up to maxParallel threads at once
    m_decompressNode.reset(new function_node<Message Message>(*m_graph, m_maxParallel, f));

    // The sequencer node enforces the correct ordering based upon the message id
    m_ordering.reset(new sequencer_node<Message>(*m_graph, [] (const Message& msg) -> unsigned int {
      return msg.id;
      }));

    // Sampling of the column to the async algorithm class
    auto g = [this] (Message msg) -> continue_msg {
      const auto betas = m_bayes->processColumnAsync(msg.marker.get());
      msg.old_beta = std::get<0>(betas);
      msg.beta = std::get<1>(betas);
      msg.deltaEps= std::get<2>(betas);//we update the deltaEps of the message
      return msg;
      };

    // Sample in parallel but with a variable maxParallel2
    m_asyncSamplingNode.reset(new function_node<Message, Message>(*m_graph, maxParallel2, g));

    // Decide whether to continue calculations or discard
    auto h = [] (decision_node::input_type input,
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

    m_decisionNode.reset(new decision_node(*m_graph, maxParallel2, h));


    // Do global computation
    auto i = [this] (Message msg) -> continue_msg {
      m_bayes->updateGlobal(msg.marker.get(), msg.old_beta, msg.beta, msg.deltaEps);
      return continue_msg();
      };
    // Use the serial policy
    m_globalComputeNode.reset(new function_node<Message>(*m_graph, serial, i));

    // Limit the number of sampler
    m_limit.reset(new limiter_node<Message>(*m_graph, serial));


    // Set up the graph topology:
    //
    // orderingNode -> decompressionNode (parallel)
    //                            |
    //                            |
    //  limitNode (serial) -> samplingNode (parallel)
    //        ^                   |
    //        |                   |
    //        |               decisionNode (parallel)
    //        |                   |
    //        |                   | keep
    //        |                   |
    //        |______________globalCompute (serial)
    //
    // Run the decompressionAndSampling node in the correct order, but do not
    // wait for the most up-to-date data.

    make_edge(*m_ordering, *m_decompressNode);
    make_edge(*m_decompressNode, *m_asyncSamplingNode);
    make_edge(*m_asyncSamplingNode, *m_decisionNode);
    // Feedback that we can now sample more columns, OR
    make_edge(output_port<0>(*m_decisionNode), m_limit->decrement);
    // pass to the global computation
    make_edge(output_port<1>(*m_decisionNode), *m_globalComputeNode);

    // Feedback that we can now sample more columns
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
