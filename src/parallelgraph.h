#ifndef DENSEPARALLELGRAPH_H
#define DENSEPARALLELGRAPH_H

#include "analysisgraph.hpp"
#include "common.h"

#include "tbb/flow_graph.h"
#include <functional>
#include <memory>
#include <Eigen/Eigen>

class BayesRBase;

struct AsyncResult;

using DecompressionToken = size_t;
using AnalysisToken = size_t;

using namespace tbb::flow;

class ParallelGraph : public AnalysisGraph
{
public:
    explicit ParallelGraph(size_t decompressionTokens, size_t analysisTokens, bool useMarkerCache);
    ~ParallelGraph();

    bool isAsynchronous() const override { return true; }

    void exec(Analysis* analysis,
              unsigned int numKeptInds,
              unsigned int numIncdSnps,
              const std::vector<unsigned int> &markerIndices) override;

    // The maximum number of decompression_node bodies which can run concurrently
    size_t decompressionNodeConcurrency() const;
    void setDecompressionNodeConcurrency(size_t c);

    // The maximum number of decompressed messages in memory
    size_t decompressionTokens() const;
    void setDecompressionTokens(size_t t);

    // The maximum number of analysis_node bodies which can run concurrently
    size_t analysisNodeConcurrency() const;
    void setAnalysisNodeConcurrency(size_t c);

    // The number of asynchronous analyses to do before forcing synchronisation
    size_t analysisTokens() const;
    void setAnalysisTokens(size_t t);

private:
    struct Message {
        unsigned int id = 0;
        unsigned int snp = 0;
        unsigned int numInds = 0;
        KernelPtr kernel = nullptr;
        ConstAsyncResultPtr result = nullptr;
    };

    using DecompressionTuple = tbb::flow::tuple<DecompressionToken, Message>;
    using AnalysisTuple = tbb::flow::tuple<AnalysisToken, DecompressionTuple>;

    std::unique_ptr<graph> m_graph;

    using decompression_join_node = join_node<DecompressionTuple, queueing>;
    std::unique_ptr<decompression_join_node> m_decompressionJoinNode;

    using decompression_node = function_node<DecompressionTuple, DecompressionTuple>;
    std::unique_ptr<decompression_node> m_decompressionNode;

    using cache_reader_node = function_node<DecompressionTuple, DecompressionTuple, lightweight>;
    std::unique_ptr<cache_reader_node> m_cacheReaderNode;

    using analysis_join_node = join_node<AnalysisTuple, queueing>;
    std::unique_ptr<analysis_join_node> m_analysisJoinNode;

    using analysis_node = function_node<AnalysisTuple, AnalysisTuple>;
    std::unique_ptr<analysis_node> m_analysisNode;

    using thread_safe_update_node = function_node<AnalysisTuple, AnalysisTuple, lightweight>;
    std::unique_ptr<thread_safe_update_node> m_threadSafeUpdateNode;

    using decision_node = multifunction_node<AnalysisTuple, tbb::flow::tuple<DecompressionToken, AnalysisToken, AnalysisTuple>>;
    std::unique_ptr<decision_node> m_decisionNode;

    using global_update_node = multifunction_node<AnalysisTuple, tbb::flow::tuple<DecompressionToken, AnalysisToken>>;
    std::unique_ptr<global_update_node> m_globalUpdateNode;

    using analysis_control_node = function_node<AnalysisToken, continue_msg, lightweight>;
    std::unique_ptr<analysis_control_node> m_analysisControlNode;

    size_t m_decompressionNodeConcurrency = tbb::flow::unlimited;
    size_t m_decompressionTokens = 40;

    size_t m_analysisNodeConcurrency = tbb::flow::unlimited;
    size_t m_analysisTokens = 20;
    size_t m_analysisTokenCount = 0;

    void queueDecompressionTokens();
    void queueAnalysisTokens();
};

#endif // DENSEPARALLELGRAPH_H
