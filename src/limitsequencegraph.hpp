#ifndef LIMITSEQUENCEGRAPH_H
#define LIMITSEQUENCEGRAPH_H

#include "analysisgraph.hpp"
#include "common.h"

#include "tbb/flow_graph.h"
#include <functional>
#include <memory>

using namespace tbb::flow;


class LimitSequenceGraph : public AnalysisGraph
{
public:
    explicit LimitSequenceGraph(size_t maxParallel = 12);
    ~LimitSequenceGraph();

    bool isAsynchronous() const override { return false; }

    void exec(Analysis* analysis,
              unsigned int numKeptInds,
              unsigned int numIncdSnps,
              const std::vector<unsigned int> &markerIndices) override;

private:
    struct Message {
        unsigned int id = 0;
        unsigned int snp = 0;
        unsigned int numInds = 0;
        KernelPtr kernel = nullptr;
    };

    std::unique_ptr<graph> m_graph;
    std::unique_ptr<function_node<Message, Message>> m_decompressNode;
    std::unique_ptr<limiter_node<Message>> m_limit;
    std::unique_ptr<sequencer_node<Message>> m_ordering;
    std::unique_ptr<sequencer_node<Message>> m_ordering2;
    std::unique_ptr<function_node<Message>> m_samplingNode;
};

#endif // LIMITSEQUENCEGRAPH_H
