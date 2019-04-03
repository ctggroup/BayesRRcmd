#ifndef PREPROCESSGRAPH_H
#define PREPROCESSGRAPH_H

#include "compression.h"
#include "data.hpp"
#include "densemarker.h"
#include "markerbuilder.h"

#include "tbb/flow_graph.h"
#include <Eigen/Eigen>

#include <array>
#include <memory>

using namespace tbb::flow;

class Data;
class IndInfo;

template <typename MarkerType>
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

        using MarkerPtr = std::shared_ptr<MarkerType>;
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

template <typename MarkerType>
PreprocessGraph<MarkerType>::PreprocessGraph(size_t maxParallel)
    : m_maxParallel(maxParallel)
    , m_graph(new graph)
{
    auto processAndCompress = [] (Message msg) -> Message {
        const auto columnSize = (msg.data->numInds + 3) >> 2;

        ifstream inStream(msg.bedFile.c_str(), ios::binary);
        if (!inStream) {
            cerr << "Error: can not open the file [" + msg.bedFile + "] to read." << endl;
            return msg;
        }

        const auto offset = 3 + (msg.startSnp * columnSize);
        inStream.seekg(offset);

        MarkerBuilder<MarkerType> builder;

        for (size_t j = msg.startSnp, chunk = 0; j < msg.data->numSnps && chunk < msg.chunkSize; ++j, ++chunk ) {
            SnpInfo *snpInfo = msg.data->snpInfoVec[j];

            if (!snpInfo->included) {
                inStream.ignore(columnSize);
                continue;
            }

            builder.initialise(j, static_cast<double>(msg.data->numInds));

            for (unsigned int i = 0; i < msg.data->numInds;) {
                char ch;
                inStream.read(&ch, 1);
                if (!inStream) {
                    cerr << "Error: problem with the BED file ... has the FAM/BIM file been changed?" << endl;
                    // Abort here?
                    return msg;
                }

                bitset<8> b = ch;
                unsigned int k = 0;

                while (k < 7 && i < msg.data->numInds) {
                    if (!msg.data->indInfoVec[i]->kept) {
                        k += 2;
                    } else {
                        const unsigned int allele1 = (!b[k++]);
                        const unsigned int allele2 = (!b[k++]);

                        builder.processAllele(i, allele1, allele2);
                    }
                    i++;
                }
            }

            builder.endColumn();
            msg.snpData.at(chunk).reset(builder.build());

            // Compress the data
            if (msg.compress) {
                const auto* marker = msg.snpData.at(chunk).get();
                msg.compressedSnpData.at(chunk) = compress(marker);

                // Delete the uncompressed snp data
                msg.snpData.at(chunk).reset();
            }
        }

        return msg;
    };

    m_processAndCompressNode.reset(new function_node<Message, Message>(*m_graph, m_maxParallel, processAndCompress));

    // The sequencer node enforces the correct ordering based upon the message id
    m_ordering.reset(new sequencer_node<Message>(*m_graph, [] (const Message& msg) -> unsigned int {
        return msg.id;
    }));

    m_ordering2.reset(new sequencer_node<Message>(*m_graph, [] (const Message& msg) -> unsigned int {
        return msg.id;
    }));

    // Control the number of messages flowing through the graph
    m_limit.reset(new limiter_node<Message>(*m_graph, m_maxParallel));

    auto writeToDisk = [this] (Message msg) -> Message {
        // Write out the preprocessed data
        if (!msg.compress) {
            for (auto& dataPtr : msg.snpData) {
                if (!dataPtr)
                    continue;

                write(dataPtr.get(), m_output.get());
            }
        } else {
            std::for_each(msg.compressedSnpData.cbegin(),
                          msg.compressedSnpData.cend(),
                          [&] (const CompressedMarker& compressed) {
                if (!compressed.buffer)
                    return;

                writeCompressedDataWithIndex(compressed.buffer.get(),
                                             compressed.size,
                                             *m_output,
                                             *m_indexOutput,
                                             m_position);
            });
        }
        return msg;
    };

    m_writeNode.reset(new function_node<Message>(*m_graph, serial, writeToDisk));

    // Set up the graph topology:
    //
    // orderingNode -> limitNode -> processAndCompressNode (parallel) -> orderingNode -> writeNode (sequential)
    //                      ^                                                               |
    //                      |_______________________________________________________________|
    //
    // This ensures that we feed the writeNode with compressed data in the correct order and signal back to
    // the parallel processAndCompressNode to keep it constantly fed. This should be a self-balancing graph.
    make_edge(*m_ordering, *m_limit);
    make_edge(*m_limit, *m_processAndCompressNode);
    make_edge(*m_processAndCompressNode, *m_ordering2);
    make_edge(*m_ordering2, *m_writeNode);

    // Feedback that we can now decompress another column
    make_edge(*m_writeNode, m_limit->decrement);
}

template <typename MarkerType>
void PreprocessGraph<MarkerType>::preprocessBedFile(const string &bedFile,
                                                    const string &preprocessedBedFile,
                                                    const string &preprocessedBedIndexFile,
                                                    bool compress,
                                                    const Data *data,
                                                    size_t chunkSize)
{
    // Reset the graph from the previous iteration. This resets the sequencer node current index etc.
    m_graph->reset();
    m_position = 0;

    // Verify prerequisites and BED file
    cout << "Preprocessing bed file: " << bedFile << ", Compress data = " << (compress ? "yes" : "no") << endl;
    if (!data) {
        cerr << "Error: Cannot preprocess data with invalid Data*" << endl;
        return;
    }
    if (chunkSize < 1) {
        cerr << "Error: chunkSize must be at least 1" << endl;
        return;
    }
    if (data->numSnps == 0) {
        cerr << "Error: No SNP is retained for analysis." << endl;
        return;
    }
    if (data->numInds == 0) {
        cerr << "Error: No individual is retained for analysis." << endl;
        return;
    }

    ifstream inStream(bedFile.c_str(), ios::binary);
    if (!inStream) {
        cerr << "Error: can not open the file [" + bedFile + "] to read." << endl;
        return;
    }

    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;

    char header[3];
    inStream.read(header, 3);
    if (!inStream || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " + bedFile << endl;
        return;
    }

    m_output = std::make_unique<std::ofstream>(preprocessedBedFile.c_str(), ios::binary);
    if (m_output->fail()) {
        cerr << "Error: Unable to open the preprocessed bed file [" + preprocessedBedFile + "] for writing." << endl;
        return;
    }

    if (compress) {
        m_indexOutput = std::make_unique<std::ofstream>(preprocessedBedIndexFile.c_str(), ios::binary);
        if (m_indexOutput->fail()) {
            cerr << "Error: Unable to open the preprocessed bed index file [" + preprocessedBedIndexFile + "] for writing." << endl;
            return;
        }
    }

    size_t msgId = 0;
    for (streamsize snp = 0; snp < data->numSnps; snp += chunkSize, ++msgId) {

        Message msg {
            msgId,
            snp,
            chunkSize,
            compress,
            bedFile,
            data,
            {chunkSize, nullptr}, // snpData
            {chunkSize, {nullptr, 0}}, // compressedData
        };

        m_ordering->try_put(msg);
    }

    inStream.clear();
    inStream.close();

    // Wait for the graph to complete
    m_graph->wait_for_all();

    // Clean up
    m_output.reset();
    if (compress)
        m_indexOutput.reset();

    cout << "Finished reading PLINK BED file." << endl;
}

#endif // PREPROCESSGRAPH_H
