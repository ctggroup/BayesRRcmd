#include "preprocessgraph.h"

#include "compression.h"
#include "data.hpp"

static std::streamsize kHeaderSize = 3;

PreprocessGraph::PreprocessGraph(size_t maxParallel)
    : m_maxParallel(maxParallel)
    , m_graph(new graph)
{
    auto processAndCompress = [] (Message msg) -> Message {
        const auto columnSize = (msg.data->numInds + 3) >> 2;
        const auto maxCompressedOutputSize = maxCompressedDataSize<double>(msg.data->numInds);

        ifstream inStream(msg.bedFile.c_str(), ios::binary);
        if (!inStream) {
            cerr << "Error: can not open the file [" + msg.bedFile + "] to read." << endl;
            return msg;
        }

        const auto offset = kHeaderSize + (msg.startSnp * columnSize);
        inStream.seekg(offset);

        for (size_t j = msg.startSnp, chunk = 0; j < msg.data->numSnps && chunk < msg.chunkSize; ++j, ++chunk ) {
            SnpInfo *snpInfo = msg.data->snpInfoVec[j];

            if (!snpInfo->included) {
                inStream.ignore(columnSize);
                continue;
            }

            // Make a note of which individuals have a missing genotype
            vector<Eigen::Index> missingIndices;
            unsigned int nmiss = 0;

            // Create some scratch space to preprocess the raw data
            msg.snpData.at(chunk) = std::make_shared<VectorXd>(msg.data->numInds);
            VectorXd& snpData = *msg.snpData.at(chunk);
            double sum = 0.0;

            for (unsigned int i = 0, ind = 0; i < msg.data->numInds;) {
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

                        if (allele1 == 0 && allele2 == 1) {  // missing genotype
                            // Don't store a marker value like this as it requires floating point comparisons later
                            // which are not done properly. Instead, store the index of the individual in a vector and simply
                            // iterate over the collected indices. Also means iterating over far fewer elements which may
                            // make a noticeable difference as this scales up.
                            missingIndices.push_back(ind++);
                            ++nmiss;
                        } else {
                            const auto value = allele1 + allele2;
                            snpData[ind++] = value;
                            sum += value;
                        }
                    }
                    i++;
                }
            }

            // Fill missing values with the mean
            const double mean = sum / double(msg.data->numInds - nmiss);
            if (nmiss) {
                for (const auto index : missingIndices)
                    snpData[index] = mean;
            }

            // Standardize genotypes
            snpData.array() -= snpData.mean();
            const auto sqn = snpData.squaredNorm();
            const auto sigma = 1.0 / (sqrt(sqn / (double(msg.data->numInds - 1))));
            snpData.array() *= sigma;

            // Compress the data
            if (msg.compress) {
                msg.compressedSnpData.at(chunk).reset(new unsigned char[maxCompressedOutputSize]);
                msg.sizes.at(chunk) = compressData(snpData,
                                               msg.compressedSnpData.at(chunk).get(),
                                               maxCompressedOutputSize);
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

                VectorXd& snpData = *dataPtr;
                m_output->write(reinterpret_cast<char *>(&snpData[0]), msg.data->numInds * sizeof(double));
            }
        } else {
            for (size_t i = 0; i < msg.compressedSnpData.size(); ++i) {
                auto& compressedData = msg.compressedSnpData.at(i);
                if (!compressedData)
                    continue; // We might not have a full chunk of data

                writeCompressedDataWithIndex(compressedData.get(),
                                 msg.sizes.at(i),
                                *m_output,
                                *m_indexOutput,
                                m_position);
            }
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

void PreprocessGraph::preprocessBedFile(const string &bedFile,
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

    char header[kHeaderSize];
    inStream.read(header, kHeaderSize);
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
            {chunkSize, nullptr}, // compressedData
            Message::SizeList(chunkSize, 0)
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
