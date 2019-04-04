#ifndef COMMON_H
#define COMMON_H

enum DataType : unsigned int {
    None = 0,
    Dense,
    SparseEigen,
    SparseRagged
};

class MarkerBuilder;
MarkerBuilder* builderForType(const DataType type);

// An entry for the index to the compressed preprocessed bed file
struct IndexEntry {
    long pos;
    long size;
};

#endif // COMMON_H
