#ifndef COMPRESSION_H
#define COMPRESSION_H

#include <Eigen/Eigen>

using namespace Eigen;

struct DataAndSize {
    unsigned char *data = nullptr;
    long size = 0;
};

DataAndSize compressData(const VectorXf &snpData);

#endif // COMPRESSION_H
