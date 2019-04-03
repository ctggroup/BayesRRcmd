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

#endif // COMMON_H
