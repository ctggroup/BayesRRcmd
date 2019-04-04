#ifndef COMMON_H
#define COMMON_H

#include <string>

enum DataType : unsigned int {
    None = 0,
    Dense,
    SparseEigen,
    SparseRagged
};

class MarkerBuilder;
MarkerBuilder* builderForType(const DataType type);

std::string ppFileForType(DataType type, const std::string &bedFile);
std::string ppIndexFileForType(DataType type, const std::string &bedFile);

// An entry for the index to the compressed preprocessed bed file
struct IndexEntry {
    long pos;
    long size;
};

#endif // COMMON_H
