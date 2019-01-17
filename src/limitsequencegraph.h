#ifndef LIMITSEQUENCEGRAPH_H
#define LIMITSEQUENCEGRAPH_H

#include "tbb/flow_graph.h"
#include <functional>
#include <memory>

class BayesRRmz;

using namespace tbb::flow;

class LimitSequenceGraph
{
public:
    LimitSequenceGraph(BayesRRmz *bayes, size_t maxParallel = 12);

    void exec(unsigned int numKeptInds,
              unsigned int numIncdSnps,
              const std::vector<unsigned int> &markerIndices);

private:
    struct Message {
        unsigned int id;
        unsigned int marker;
        unsigned int numKeptInds;
        unsigned char *decompressBuffer = nullptr;
    };

    BayesRRmz *m_bayes;
    size_t m_maxParallel;
    std::unique_ptr<graph> m_graph;
    std::unique_ptr<function_node<Message, Message>> m_decompressNode;
    std::unique_ptr<limiter_node<Message>> m_limit;
    std::unique_ptr<sequencer_node<Message>> m_ordering;
    std::unique_ptr<sequencer_node<Message>> m_ordering2;
    std::unique_ptr<function_node<Message>> m_samplingNode;
};

#endif // LIMITSEQUENCEGRAPH_H
