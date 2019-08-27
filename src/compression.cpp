#include "compression.h"

#include <iostream>

int prepareStream(z_stream &strm) {
    strm.zalloc = nullptr;
    strm.zfree = nullptr;
    strm.opaque = nullptr;
    const int level = -1;
    return deflateInit(&strm, level);
}

unsigned long compressData(z_stream &strm, unsigned long outputSize) {
    const int flush = Z_FINISH;
    if (deflate(&strm, flush) != Z_STREAM_END) {
        std::cout << "Error compressing data" << std::endl;
        return 0;
    }
    const auto compressedSize = outputSize - strm.avail_out;
    // std::cout << "compressedSize = " << compressedSize << " bytes = "
    //           << compressedSize / 1024 << " KiB" << std::endl;

    // Clean up
    (void) deflateEnd(&strm);

    // DEBUG: Verify compressed data can be decompressed to reproduce the original data
    /*
    z_stream strm2;
    strm2.zalloc = nullptr;
    strm2.zfree = nullptr;
    strm2.opaque = nullptr;
    strm2.avail_in = 0;
    strm2.next_in = nullptr;
    ret = inflateInit(&strm2);
    if (ret != Z_OK) {
        std::cout << "Failed to verify compressed data" << std::endl;
        return compressedSize;
    }
    unsigned char *checkBuffer = new unsigned char[inputSize];
    strm2.next_out = checkBuffer;
    strm2.avail_out = inputSize;
    strm2.next_in = outputBuffer;
    strm2.avail_in = static_cast<unsigned int>(compressedSize);
    ret = inflate(&strm2, flush);
    if (ret != Z_STREAM_END) {
        std::cout << "Failed to verify compressed data" << std::endl;
        return compressedSize;
    }
    // Compare input and re-extracted data
    {
        Map<VectorXf> decompressedSnpData(reinterpret_cast<float *>(checkBuffer), snpData.size());
        for (int i = 0; i < snpData.size(); ++i) {
            const auto delta = snpData[i] - decompressedSnpData[i];
            std::cout << i << ": delta = " << delta << std::endl;
        }
    }
    // Cleanup
    delete[] checkBuffer;
    (void) inflateEnd(&strm2);
    */

    return compressedSize;
}

unsigned long compressData(const VectorXd &snpData, unsigned char *outputBuffer, unsigned long outputSize)
{
    // Initialise zlib
    z_stream strm;
    if (prepareStream(strm) != Z_OK)
        return 0;

    // Compress the data
    const unsigned int inputSize = static_cast<unsigned int>(snpData.size()) * sizeof(double);
    strm.avail_in = inputSize;
    strm.next_in = reinterpret_cast<unsigned char *>(const_cast<double*>(&snpData[0]));
    strm.avail_out = static_cast<unsigned int>(outputSize);
    strm.next_out = outputBuffer;

    return compressData(strm, outputSize);
}

unsigned long compressData(char *inputBuffer,
                           unsigned int inputSize,
                           unsigned char *outputBuffer,
                           unsigned long outputSize)
{
    // Initialise zlib
    z_stream strm;
    if (prepareStream(strm) != Z_OK)
        return 0;

    // Compress the data
    strm.avail_in = inputSize;
    strm.next_in = reinterpret_cast<unsigned char *>(inputBuffer);
    strm.avail_out = static_cast<unsigned int>(outputSize);
    strm.next_out = outputBuffer;

    return compressData(strm, outputSize);
}

void extractData(unsigned char *compressedData,
                 unsigned int compressedDataSize,
                 unsigned char *outputBuffer,
                 unsigned int outputBufferSize)
{
    z_stream strm;
    strm.zalloc = nullptr;
    strm.zfree = nullptr;
    strm.opaque = nullptr;
    strm.avail_in = 0;
    strm.next_in = nullptr;
    auto ret = inflateInit(&strm);
    if (ret != Z_OK)
        throw("Failed to verify compressed data");

    strm.next_out = outputBuffer;
    strm.avail_out = outputBufferSize;
    strm.next_in = compressedData;
    strm.avail_in = compressedDataSize;
    const int flush = Z_FINISH;
    ret = inflate(&strm, flush);
    if (ret != Z_STREAM_END)
        throw("Failed to verify compressed data");

    (void) inflateEnd(&strm);
}

void writeUncompressedDataWithIndex(const unsigned char *data,
                                    const unsigned long size,
                                    std::ostream &outStream,
                                    std::ostream &indexStream,
                                    unsigned long &pos)
{
    writeCompressedDataWithIndex(data, size, size, outStream, indexStream, pos);
}

void writeCompressedDataWithIndex(const unsigned char *data,
                                  const unsigned long compressedSize,
                                  const unsigned long originalSize,
                                  std::ostream &outStream,
                                  std::ostream &indexStream,
                                  unsigned long &pos)
{
    outStream.write(reinterpret_cast<const char *>(data),
                    static_cast<std::streamsize>(compressedSize));

    indexStream.write(reinterpret_cast<char *>(&pos),
                      sizeof(unsigned long));
    indexStream.write(reinterpret_cast<const char *>(&compressedSize),
                      sizeof(unsigned long));
    indexStream.write(reinterpret_cast<const char *>(&originalSize),
                      sizeof(unsigned long));
    pos += compressedSize;
}

void compressAndWriteWithIndex(const VectorXd &data,
                               std::ostream &outStream,
                               std::ostream &indexStream,
                               unsigned long &pos,
                               unsigned char *compressedBuffer,
                               const unsigned long maxCompressedOutputSize)
{
    const unsigned long compressedSize = compressData(data,
                                                      compressedBuffer,
                                                      maxCompressedOutputSize);

    const unsigned long originalSize = data.size() * sizeof(double);
    writeCompressedDataWithIndex(compressedBuffer,
                                 compressedSize,
                                 originalSize,
                                 outStream,
                                 indexStream,
                                 pos);
}
