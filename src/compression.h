#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <Eigen/Eigen>
#include <zlib.h>

using namespace Eigen;

template <typename T>
unsigned long maxCompressedDataSize(const unsigned int count)
{
    // Initialise zlib
    z_stream strm;
    strm.zalloc = nullptr;
    strm.zfree = nullptr;
    strm.opaque = nullptr;
    const int level = -1;
    auto ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return 0;

    // Calculate the maximum buffer size needed to hold the compressed data
    const unsigned int inputSize = count * sizeof(T);
    strm.avail_in = inputSize;
    const auto maxOutputSize = deflateBound(&strm, inputSize);
    //std::cout << "maxSize = " << maxOutputSize << " bytes = " << maxOutputSize / 1024 << " KiB" << std::endl;

    // Clean up
    (void) deflateEnd(&strm);

    return maxOutputSize;
}

int prepareStream(z_stream &strm);
unsigned long compressData(z_stream &strm, unsigned long outputSize);

unsigned long compressData(const VectorXd &snpData,
                           unsigned char *outputBuffer,
                           unsigned long outputSize);

unsigned long compressData(char *inputBuffer,
                           unsigned int inputSize,
                           unsigned char *outputBuffer,
                           unsigned long outputSize);

template<typename T>
unsigned long compressData(const std::vector<T> &vector,
                           unsigned char *outputBuffer,
                           unsigned long outputSize)
{
    // Initialise zlib
    z_stream strm;
    if (prepareStream(strm) != Z_OK)
        return 0;

    // Compress the data
    const unsigned int inputSize = static_cast<unsigned int>(vector.size()) * sizeof(T);
    strm.avail_in = inputSize;
    strm.next_in = reinterpret_cast<unsigned char *>(const_cast<T*>(vector.data()));
    strm.avail_out = static_cast<unsigned int>(outputSize);
    strm.next_out = outputBuffer;

    return compressData(strm, outputSize);
}

void writeUncompressedDataWithIndex(const unsigned char *data,
                                    const unsigned long size,
                                    std::ostream &outStream,
                                    std::ostream &indexStream,
                                    unsigned long &pos);

void writeCompressedDataWithIndex(const unsigned char *data,
                                  const unsigned long compressedSize,
                                  const unsigned long originalSize,
                                  std::ostream &outStream,
                                  std::ostream &indexStream,
                                  unsigned long &pos);

void compressAndWriteWithIndex(const VectorXd &data,
        std::ostream &outStream,
        std::ostream &indexStream,
        unsigned long &pos,
        unsigned char *compressedBuffer,
        const unsigned long maxCompressedOutputSize);

void extractData(unsigned char *compressedData,
                 unsigned int compressedDataSize,
                 unsigned char *outputBuffer,
                 unsigned int outputBufferSize);

#endif // COMPRESSION_H
