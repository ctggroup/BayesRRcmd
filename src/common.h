#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <iostream>

enum class InputType : unsigned int {
    Unknown = 0,
    BED,
    CSV
};

enum class PreprocessDataType : unsigned int {
    None = 0,
    Dense,
    SparseEigen,
    SparseRagged
};

std::ostream &operator<<(std::ostream &os, const PreprocessDataType &obj);

class MarkerBuilder;
MarkerBuilder* builderForType(const PreprocessDataType type);

std::string fileWithSuffix(const std::string &dataFile, const std::string &suffix);

std::string ppFileForType(PreprocessDataType type, const std::string &dataFile);
std::string ppIndexFileForType(PreprocessDataType type, const std::string &dataFile);

InputType getInputType(const std::string &dataFile);

// An entry for the index to the compressed preprocessed bed file
struct IndexEntry {
    unsigned long pos = 0;
    unsigned long compressedSize = 0;
    unsigned long originalSize = 0;
};

//template to read csv files into an eigen vector.
template<typename M>
M load_csv (const std::string & path);


#endif // COMMON_H
