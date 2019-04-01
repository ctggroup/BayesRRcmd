#ifndef PREPROCESSGRAPH_H
#define PREPROCESSGRAPH_H

#include "tbb/flow_graph.h"
#include <Eigen/Eigen>

#include <array>
#include <memory>

using namespace tbb::flow;

class Data;
class IndInfo;

class PreprocessGraph
{
public:
    PreprocessGraph(size_t maxParallel);

    void preprocessBedFile(const std::string &bedFile, const std::string &preprocessedBedFile, const std::string &preprocessedBedIndexFile, bool compress, const Data *data, size_t chunkSize);

protected:
    struct Message {
        size_t id = 0;
        std::streamsize startSnp = 0;
        size_t chunkSize = 0;
        bool compress = false;

        std::string bedFile = "";

        const Data* data = nullptr;

        using SnpData = Eigen::VectorXd;
        using SnpDataPtr = std::shared_ptr<SnpData>;
        using SnpDataPtrList = std::vector<SnpDataPtr>;
        SnpDataPtrList snpData;

        using CompressedData = unsigned char[];
        using CompressedDataPtr = std::shared_ptr<CompressedData>;
        using CompressedDataPtrList = std::vector<CompressedDataPtr>;
        CompressedDataPtrList compressedSnpData;

        using SizeList = std::vector<unsigned long>;
        SizeList sizes;
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
