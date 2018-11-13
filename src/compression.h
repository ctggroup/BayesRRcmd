#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <Eigen/Eigen>

using namespace Eigen;

unsigned long maxCompressedDataSize(const unsigned int numFloats);

unsigned long compressData(const VectorXf &snpData,
                           unsigned char *outputBuffer,
                           unsigned long outputSize);

#endif // COMPRESSION_H
