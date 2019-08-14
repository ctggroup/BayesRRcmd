#ifndef PREPROCESSGRAPH_H
#define PREPROCESSGRAPH_H

#include "common.h"
#include "compression.h"
#include "data.hpp"
#include "marker.h"

#include "tbb/flow_graph.h"
#include <Eigen/Eigen>

#include <array>
#include <memory>

using namespace tbb::flow;

class PreprocessGraph
{
public:
    explicit PreprocessGraph(size_t maxParallel);
    ~PreprocessGraph();

    void preprocessBedFile(const std::string &dataFile,
                           const PreprocessDataType type,
                           const bool compress,
                           const Data *data,
                           const size_t chunkSize);

protected:
    struct Message {
        PreprocessDataType type = PreprocessDataType::None;
        size_t id = 0;
        std::streamsize startSnp = 0;
        size_t chunkSize = 0;
        bool compress = false;

        std::string bedFile = "";

        const Data* data = nullptr;

        using MarkerPtrList = std::vector<MarkerPtr>;
        MarkerPtrList snpData;

        using CompressedMarkerList = std::vector<CompressedMarker>;
        CompressedMarkerList compressedSnpData;
    };

    size_t m_maxParallel = 1;
    std::unique_ptr<graph> m_graph;
    std::unique_ptr<function_node<Message, Message>> m_processAndCompressNode;
    std::unique_ptr<limiter_node<Message>> m_limit;
    std::unique_ptr<sequencer_node<Message>> m_ordering;
    std::unique_ptr<sequencer_node<Message>> m_ordering2;
    std::unique_ptr<function_node<Message>> m_writeNode;

    std::unique_ptr<std::ofstream> m_output = nullptr;
    std::unique_ptr<std::ofstream> m_indexOutput = nullptr;
    unsigned long m_position = 0;
};

#endif // PREPROCESSGRAPH_H
