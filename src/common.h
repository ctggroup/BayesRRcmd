#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

enum class AnalysisType : unsigned int {
    Unknown = 0,
    Preprocess,
    PpBayes,
    AsyncPpBayes,
    Gauss,
    AsyncGauss,
    Split
};

std::ostream &operator<<(std::ostream &os, const AnalysisType &obj);

enum class InputType : unsigned int {
    Unknown = 0,
    BED,
    CSV
};

std::ostream &operator<<(std::ostream &os, const InputType &obj);

enum class PreprocessDataType : unsigned int {
    None = 0,
    Dense,
    SparseEigen,
    SparseRagged
};

std::ostream &operator<<(std::ostream &os, const PreprocessDataType &obj);

struct Marker;
using MarkerPtr = std::shared_ptr<Marker>;
using ConstMarkerPtr = std::shared_ptr<const Marker>;

struct Kernel;
using KernelPtr = std::shared_ptr<Kernel>;

struct AsyncResult;
using AsyncResultPtr = std::shared_ptr<AsyncResult>;
using ConstAsyncResultPtr = std::shared_ptr<const AsyncResult>;

class MarkerBuilder;
MarkerBuilder* builderForType(const PreprocessDataType type);

std::string fileWithSuffix(const std::string &dataFile, const std::string &suffix);

InputType getInputType(const std::string &dataFile);

// An entry for the index to the compressed preprocessed bed file
struct IndexEntry {
    unsigned long pos = 0;
    unsigned long compressedSize = 0;
    unsigned long originalSize = 0;

    auto asTuple() const { return std::tie(pos, compressedSize, originalSize); }
};

bool operator==(const IndexEntry &lhs, const IndexEntry &rhs);
bool operator!=(const IndexEntry &lhs, const IndexEntry &rhs);

using MarkerIndexList = std::vector<unsigned int>;

struct MarkerSubset;

//template to read csv files into an eigen vector.
template<typename M>
M load_csv (const std::string & path);


#endif // COMMON_H
