#include "common.h"

#include "densemarkerbuilder.h"
#include "eigensparsemarkerbuilder.h"
#include "raggedsparsemarkerbuilder.h"

#include <cassert>
#include <iostream>

#include <Eigen/Eigen>
#include <vector>
#include <fstream>
#include <type_traits>

using namespace Eigen;

std::ostream &operator<<(std::ostream &os, const AnalysisType &obj)
{
    os << static_cast<std::underlying_type<AnalysisType>::type>(obj);
    return os;
}

std::ostream &operator<<(std::ostream &os, const InputType &obj)
{
    os << static_cast<std::underlying_type<InputType>::type>(obj);
    return os;
}

std::ostream &operator<<(std::ostream &os, const PreprocessDataType &obj)
{
   os << static_cast<std::underlying_type<PreprocessDataType>::type>(obj);
   return os;
}

MarkerBuilder *builderForType(const PreprocessDataType type)
{
    switch (type)
    {
    case PreprocessDataType::Dense:
        return new DenseMarkerBuilder;

    case PreprocessDataType::SparseEigen:
        return new EigenSparseMarkerBuilder;

    case PreprocessDataType::SparseRagged:
        return new RaggedSparseMarkerBuilder;

    case PreprocessDataType::None:
        // Fall through
    default:
        std::cerr << "builderForType - unhandled DataType:" << type << std::endl;
        assert(false);
        return nullptr;
    }
}

std::string fileWithoutExtension(const std::string &dataFile)
{
    auto const pos = dataFile.find_last_of('.');
    if (pos == std::string::npos) {
        std::cout << "Could remove file extension: " << dataFile << std::endl;
        return dataFile;
    }
    return dataFile.substr(0, pos);
}

std::string fileWithSuffix(const std::string &dataFile, const std::string &suffix)
{
    return fileWithoutExtension(dataFile) + suffix;
}

std::string ppFileForType(PreprocessDataType type, const std::string &dataFile)
{
    const auto fileName = fileWithoutExtension(dataFile);

    switch (type) {
    case PreprocessDataType::Dense:
        return fileName + ".ppbed";

    case PreprocessDataType::SparseEigen:
        return fileName +  ".eigen.sparsebed";

    case PreprocessDataType::SparseRagged:
        return fileName +  ".ragged.sparsebed";

    default:
        std::cerr << "ppFileForType - unsupported DataType: "
             << type
             << std::endl;
        assert(false);
        return {};
    }
}

std::string ppIndexFileForType(PreprocessDataType type, const std::string &dataFile)
{
    const auto fileName = fileWithoutExtension(dataFile);

    switch (type) {
    case PreprocessDataType::Dense:
        return fileName +  ".ppbedindex";

    case PreprocessDataType::SparseEigen:
        return fileName +  ".eigen.sparsebedindex";

    case PreprocessDataType::SparseRagged:
        return fileName + ".ragged.sparsebedindex";

    default:
        std::cerr << "ppIndexFileForType - unsupported DataType: "
             << type
             << std::endl;
        assert(false);
        return {};
    }
}

InputType getInputType(const std::string &dataFile)
{
    auto endsWith = [](std::string const & value, std::string const & ending)
    {
        if (ending.size() > value.size())
            return false;

        return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    };

    if (endsWith(dataFile, ".bed"))
        return InputType::BED;
    else if (endsWith(dataFile, ".csv"))
        return InputType::CSV;
    else
        return InputType::Unknown;
}



template<typename M>
M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), rows, values.size()/rows);
}
