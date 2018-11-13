#include "compression.h"

#include <zlib.h>
#include <iostream>

DataAndSize compressData(const VectorXf &snpData)
{
    DataAndSize result;

    // Initialise zlib
    z_stream strm;
    strm.zalloc = nullptr;
    strm.zfree = nullptr;
    strm.opaque = nullptr;
    const int level = -1;
    auto ret = deflateInit(&strm, level);
    if (ret != Z_OK)
        return result;

    // Calculate the maximum buffer size needed to hold the compressed data
    const unsigned int inputSize = static_cast<unsigned int>(snpData.size()) * sizeof(float);
    strm.avail_in = inputSize;
    strm.next_in = reinterpret_cast<unsigned char *>(const_cast<float*>(&snpData[0]));
    const auto maxOutputSize = deflateBound(&strm, inputSize);
    //std::cout << "maxSize = " << maxOutputSize << " bytes = " << maxOutputSize / 1024 << " KiB" << std::endl;

    // Create a suitable output buffer
    unsigned char *output = new unsigned char[maxOutputSize];
    strm.avail_out = static_cast<unsigned int>(maxOutputSize);
    strm.next_out = output;

    // Compress the data
    const int flush = Z_FINISH;
    ret = deflate(&strm, flush);
    if (ret != Z_STREAM_END) {
        std::cout << "Error compressing data" << std::endl;
        return result;
    }
    const auto compressedSize = maxOutputSize - strm.avail_out;
    std::cout << "compressedSize = " << compressedSize << " bytes = "
              << compressedSize / 1024 << " KiB" << std::endl;

    // Clean up
    (void) deflateEnd(&strm);

    result.data = output;
    result.size = long(compressedSize);
    return result;
}
