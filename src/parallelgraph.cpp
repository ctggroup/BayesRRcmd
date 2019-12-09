#include "parallelgraph.h"

#include "compression.h"
#include "BayesRBase.hpp"
#include "hybridmpisyncmanager.h"
#include "kernel.h"
#include "markerbuilder.h"
#include "markercache.h"

#include <iostream>

ParallelGraph::ParallelGraph(size_t maxDecompressionTokens, size_t maxAnalysisTokens, bool useMarkerCache, bool useMpi)
    : AnalysisGraph()
    , m_graph(new graph)
    , m_mpiGraph(new graph)
    , m_syncManager(new HybridMpi::SyncManager(useMpi))
    , m_decompressionTokens(maxDecompressionTokens)
    , m_analysisTokens(maxAnalysisTokens)
{
    m_decompressionJoinNode.reset(new decompression_join_node(*m_graph));
    m_mpiDecompressionJoinNode.reset(new decompression_join_node(*m_mpiGraph));

    auto diskReader = [this] (DecompressionTuple tuple) -> DecompressionTuple {
        auto &msg = std::get<1>(tuple);
        // Read the column from disk
        std::unique_ptr<MarkerBuilder> builder{m_analysis->markerBuilder()};
        builder->initialise(msg.snp, m_numInds);
        const auto index = m_analysis->indexEntry(msg.snp);
        if (m_analysis->compressed()) {
            builder->decompress(m_analysis->compressedData(), index);
        } else {
            builder->read(m_analysis->preprocessedFile(), index);
        }
        msg.kernel = m_analysis->kernelForMarker(builder->build());
        return tuple;
    };

    if (useMarkerCache) {
        auto cacheReader = [this] (DecompressionTuple tuple) -> DecompressionTuple {
            auto &msg = std::get<1>(tuple);
            msg.kernel = m_analysis->kernelForMarker(markerCache()->marker(msg.snp));
            return tuple;
        };

        m_cacheReaderNode.reset(new cache_reader_node(*m_graph,
                                                      m_decompressionNodeConcurrency,
                                                      cacheReader));
    } else {
        m_decompressionNode.reset(new decompression_node(*m_graph,
                                                         m_decompressionNodeConcurrency,
                                                         diskReader));
    }

    m_mpiDecompressionNode.reset(new decompression_node(*m_mpiGraph,
                                                        m_decompressionNodeConcurrency,
                                                        diskReader));

    m_analysisJoinNode.reset(new analysis_join_node(*m_graph));
    m_mpiAnalysisJoinNode.reset(new analysis_join_node(*m_mpiGraph));

    // Sampling of the column to the async algorithm class
    auto processColumn = [this] (AnalysisTuple tuple) -> AnalysisTuple {
        auto &msg = std::get<1>(std::get<1>(tuple));
        msg.result = m_analysis->processColumnAsync(msg.kernel);
        return tuple;
    };

    m_analysisNode.reset(new analysis_node(*m_graph, m_analysisNodeConcurrency, processColumn));

    auto processResult = [this] (AnalysisTuple tuple) -> AnalysisTuple {
        auto &msg = std::get<1>(std::get<1>(tuple));
        const auto result = m_syncManager->resultForMarker(msg.snp);
        m_analysis->processResult(msg.kernel, result);
        msg.result = result;
        return tuple;
    };

    m_mpiAnalysisNode.reset(new analysis_node(*m_mpiGraph, m_analysisNodeConcurrency, processResult));

    auto threadSafeUpdate = [this] (AnalysisTuple tuple) -> AnalysisTuple {
        auto &msg = std::get<1>(std::get<1>(tuple));
        m_analysis->doThreadSafeUpdates(msg.result);
        return tuple;
    };

    m_threadSafeUpdateNode.reset(new thread_safe_update_node(*m_graph,
                                                             serial,
                                                             threadSafeUpdate));

    // Decide whether to continue calculations or discard
    auto h = [] (decision_node::input_type input,
            decision_node::output_ports_type &outputPorts) {

        auto &decompressionTuple = std::get<1>(input);
        auto &msg = std::get<1>(decompressionTuple);

        if (msg.result->betaOld != 0.0 || msg.result->beta != 0.0) {
            // Do global computation
            std::get<2>(outputPorts).try_put(input);
        } else {
            // Discard
            std::get<0>(outputPorts).try_put(std::get<0>(decompressionTuple));
            std::get<1>(outputPorts).try_put(std::get<0>(input));
        }
    };

    m_decisionNode.reset(new decision_node(*m_graph, unlimited, h));

    // Do global computation
    auto globalUpdate = [this, useMpi] (global_update_node::input_type input,
            global_update_node::output_ports_type &outputPorts) {

        auto &decompressionTuple = std::get<1>(input);
        auto &msg = std::get<1>(decompressionTuple);

        m_analysis->updateGlobal(msg.kernel, msg.result);

        if (useMpi)
            m_syncManager->accumulate(msg.kernel, msg.result);

        std::get<0>(outputPorts).try_put(std::get<0>(decompressionTuple));
        std::get<1>(outputPorts).try_put(std::get<0>(input));
    };
    // Use the serial policy
    m_globalUpdateNode.reset(new global_update_node(*m_graph, serial, globalUpdate));

    auto mpiGlobalUpdate = [this] (global_update_node::input_type input,
            global_update_node::output_ports_type &outputPorts) {

        auto &decompressionTuple = std::get<1>(input);
        auto &msg = std::get<1>(decompressionTuple);

        m_analysis->updateGlobal(msg.kernel, msg.result);

        std::get<0>(outputPorts).try_put(std::get<0>(decompressionTuple));
        std::get<1>(outputPorts).try_put(std::get<0>(input));
    };
    // Use the serial policy
    m_mpiGlobalUpdateNode.reset(new global_update_node(*m_mpiGraph, serial, mpiGlobalUpdate));

    if (useMpi) {
        // Force synchronisation after m_anaylsisToken analyses
        auto doUpdateMpi = [this](AnalysisToken t) -> continue_msg {
            (void) t; // Unused
            --m_analysisTokenCount;

            if (m_analysisTokenCount == 0) {
                // Synchronise values from other processes
                const auto markers = m_syncManager->synchroniseMarkers();
                if (!markers.empty()) {
                    // Populate the analysis join node with enough tokens to work through our markers
                    // for(AnalysisToken t = 0; t < markers.size(); ++t)
                    //     input_port<0>(*m_mpiAnalysisJoinNode).try_put(t);

                    // Load from disk and calculate epsilon updates
                    for (unsigned int i = 0; i < markers.size(); ++i) {
                        // Populate the analysis join node with a token for each marker
                        input_port<0>(*m_mpiAnalysisJoinNode).try_put(i);

                        Message msg = { i, markers.at(i), nullptr, nullptr };
                        input_port<1>(*m_mpiDecompressionJoinNode).try_put(msg);
                    }

                    // Wait for the graph to complete
                    m_mpiGraph->wait_for_all();
                }
                m_syncManager->reset();

                // Allow the next set of analyses to take place
                queueAnalysisTokens();
            }

            return continue_msg();
        };
        m_mpiAnalysisControlNode.reset(new mpi_analysis_control_node(*m_graph, serial, doUpdateMpi));
    } else {
        // Force synchronisation after m_anaylsisToken analyses
        auto j = [this](AnalysisToken t) -> continue_msg {
            (void) t; // Unused
            --m_analysisTokenCount;

            if (m_analysisTokenCount == 0) {
                // Allow the next set of analyses to take place
                queueAnalysisTokens();
            }

            return continue_msg();
        };
        m_analysisControlNode.reset(new analysis_control_node(*m_graph, serial, j));
    }

    m_mpiInternalAnalysisControlNode.reset(new analysis_control_node(*m_mpiGraph, serial, [](AnalysisToken t) -> continue_msg {
        return continue_msg();
    }));

    // Set up the graph topology:
#if defined(TBB_PREVIEW_FLOW_GRAPH_TRACE)
    m_graph->set_name("ParallelGraph");
    m_analysisNode->set_name("analysis_node");
    m_decisionNode->set_name("decision_node");
    m_globalUpdateNode->set_name("global_update_node");
    m_decompressionJoinNode->set_name("decompression_join_node");
    if(useMarkerCache)
        m_cacheReaderNode->set_name("cache_reader_node");
    else
        m_decompressionNode->set_name("decompression_node");
    m_analysisJoinNode->set_name("analysis_join_node");
    m_analysisControlNode->set_name("analysis_control_node");
#endif

    // Parallel graph
    if (useMarkerCache) {
        make_edge(*m_decompressionJoinNode, *m_cacheReaderNode);
        make_edge(*m_cacheReaderNode, input_port<1>(*m_analysisJoinNode));
    } else {
        make_edge(*m_decompressionJoinNode, *m_decompressionNode);
        make_edge(*m_decompressionNode, input_port<1>(*m_analysisJoinNode));
    }
    make_edge(*m_analysisJoinNode, *m_analysisNode);
    make_edge(*m_analysisNode, *m_threadSafeUpdateNode);
    make_edge(*m_threadSafeUpdateNode, *m_decisionNode);
    make_edge(output_port<0>(*m_decisionNode), input_port<0>(*m_decompressionJoinNode));
    make_edge(output_port<2>(*m_decisionNode), *m_globalUpdateNode);
    make_edge(output_port<0>(*m_globalUpdateNode), input_port<0>(*m_decompressionJoinNode));
    if (useMpi) {
        make_edge(output_port<1>(*m_decisionNode), *m_mpiAnalysisControlNode);
        make_edge(output_port<1>(*m_globalUpdateNode), *m_mpiAnalysisControlNode);
    } else {
        make_edge(output_port<1>(*m_decisionNode), *m_analysisControlNode);
        make_edge(output_port<1>(*m_globalUpdateNode), *m_analysisControlNode);
    }

    // internal MPI graph
    make_edge(*m_mpiDecompressionJoinNode, *m_mpiDecompressionNode);
    make_edge(*m_mpiDecompressionNode, input_port<1>(*m_mpiAnalysisJoinNode));
    make_edge(*m_mpiAnalysisJoinNode, *m_mpiAnalysisNode);
    make_edge(*m_mpiAnalysisNode, *m_mpiGlobalUpdateNode);
    make_edge(output_port<0>(*m_mpiGlobalUpdateNode), input_port<0>(*m_mpiDecompressionJoinNode));
    make_edge(output_port<1>(*m_mpiGlobalUpdateNode), *m_mpiInternalAnalysisControlNode);
}

ParallelGraph::~ParallelGraph()
{
    m_graph->wait_for_all();
}

void ParallelGraph::exec(Analysis *analysis,
                              unsigned int numInds,
                              const std::vector<unsigned int> &markerIndices)
{
    if (!analysis) {
        std::cerr << "Cannot run ParallelGraph without bayes" << std::endl;
        return;
    }

    // Set our Bayes for this run
    m_analysis = analysis;
    m_numInds = numInds;
    const auto numSnps = static_cast<unsigned int>(markerIndices.size());
    m_markersRemaining = numSnps;

    // Do not allow Eigen to parallalize during ParallelGraph execution.
    const auto eigenThreadCount = Eigen::nbThreads();
    Eigen::setNbThreads(0);

    // Reset the graph from the previous iteration.
    m_graph->reset();
    queueDecompressionTokens();
    queueAnalysisTokens();

    // Push some messages into the top of the graph to be processed - representing the column indices
    for (unsigned int i = 0; i < numSnps; ++i) {
        Message msg = { i, markerIndices[i], nullptr, nullptr };
        input_port<1>(*m_decompressionJoinNode).try_put(msg);
    }

    // Wait for the graph to complete
    m_graph->wait_for_all();

    // Turn Eigen threading back on.
    Eigen::setNbThreads(eigenThreadCount);

    // Clean up
    m_analysis = nullptr;
}

size_t ParallelGraph::decompressionNodeConcurrency() const
{
    return m_decompressionNodeConcurrency;
}

void ParallelGraph::setDecompressionNodeConcurrency(size_t c)
{
    m_decompressionNodeConcurrency = c;
}

size_t ParallelGraph::decompressionTokens() const
{
    return m_decompressionTokens;
}

void ParallelGraph::setDecompressionTokens(size_t t)
{
    m_decompressionTokens = t;
}

size_t ParallelGraph::analysisNodeConcurrency() const
{
    return m_analysisNodeConcurrency;
}

void ParallelGraph::setAnalysisNodeConcurrency(size_t c)
{
    m_analysisNodeConcurrency = c;
}

size_t ParallelGraph::analysisTokens() const
{
    return m_analysisTokens;
}

void ParallelGraph::setAnalysisTokens(size_t t)
{
    m_analysisTokens = t;
}

void ParallelGraph::queueDecompressionTokens()
{
    for(DecompressionToken t = 0; t < m_decompressionTokens; ++t) {
        input_port<0>(*m_decompressionJoinNode).try_put(t);
        input_port<0>(*m_mpiDecompressionJoinNode).try_put(t);
    }
}

void ParallelGraph::queueAnalysisTokens()
{
    const size_t tokens = m_markersRemaining > m_analysisTokens ? m_analysisTokens : m_markersRemaining;
    for(AnalysisToken t = 0; t < tokens; ++t)
        input_port<0>(*m_analysisJoinNode).try_put(t);

    m_analysisTokenCount = tokens;
    m_markersRemaining -= tokens;
}
