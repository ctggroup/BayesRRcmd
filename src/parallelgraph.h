#ifndef DENSEPARALLELGRAPH_H
#define DENSEPARALLELGRAPH_H

#include "analysisgraph.hpp"

#include "tbb/flow_graph.h"
#include <functional>
#include <memory>
#include <Eigen/Eigen>

class BayesRBase;

struct Marker;

using namespace tbb::flow;

class ParallelGraph : public AnalysisGraph
{
public:
    explicit ParallelGraph(size_t maxParallel = 6);

    bool isAsynchronous() const override { return true; }

    void exec(BayesRBase* bayes,
              unsigned int numKeptInds,
              unsigned int numIncdSnps,
              const std::vector<unsigned int> &markerIndices) override;

private:
    struct Message {
        unsigned int id = 0;
        unsigned int snp = 0;
        unsigned int numInds = 0;
        std::shared_ptr<Marker> marker = nullptr;
        Eigen::VectorXd deltaEps; //vector that stores the epsilon update only
        double old_beta = 0.0;
        double beta = 0.0;
    };

    std::unique_ptr<graph> m_graph;
    std::unique_ptr<function_node<Message, Message>> m_decompressNode;
    std::unique_ptr<function_node<Message, Message>> m_asyncSamplingNode;
    std::unique_ptr<limiter_node<Message>> m_limit;
    std::unique_ptr<sequencer_node<Message>> m_ordering;

    using decision_node = multifunction_node<Message, tbb::flow::tuple<continue_msg, Message> >;
    std::unique_ptr<decision_node> m_decisionNode;
    std::unique_ptr<function_node<Message>> m_globalComputeNode;
};

#endif // DENSEPARALLELGRAPH_H
