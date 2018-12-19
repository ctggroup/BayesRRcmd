#ifndef LIMITSEQUENCEGRAPH_H
#define LIMITSEQUENCEGRAPH_H

#include "tbb/flow_graph.h"
#include <memory>

using namespace tbb::flow;

class LimitSequenceGraph
{
public:
    LimitSequenceGraph(size_t maxParallel = 4);

    void exec();

private:
    struct Message {
        int id;
        float *data;
    };

    size_t m_maxParallel;
    std::unique_ptr<graph> m_graph;
    std::unique_ptr<function_node<Message, Message>> m_process;
    std::unique_ptr<limiter_node<Message>> m_limit;
    std::unique_ptr<sequencer_node<Message>> m_ordering;
    std::unique_ptr<sequencer_node<Message>> m_ordering2;
    std::unique_ptr<function_node<Message>> m_writer;
};

#endif // LIMITSEQUENCEGRAPH_H
