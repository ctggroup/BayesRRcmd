#include "common.h"

#include "densemarkerbuilder.h"
#include "eigensparsemarkerbuilder.h"
#include "raggedsparsemarkerbuilder.h"

#include <cassert>
#include <iostream>

#include <Eigen/Eigen>
#include <vector>
#include <fstream>

using namespace Eigen;

MarkerBuilder *builderForType(const DataType type)
{
    switch (type)
    {
    case DataType::Dense:
        return new DenseMarkerBuilder;

    case DataType::SparseEigen:
        return new EigenSparseMarkerBuilder;

    case DataType::SparseRagged:
        return new RaggedSparseMarkerBuilder;

    case DataType::None:
        // Fall through
    default:
        std::cerr << "builderForType - unhandled DataType:" << type << std::endl;
        assert(false);
        return nullptr;
    }
}

std::string ppFileForType(DataType type, const std::string &bedFile)
{
    switch (type) {
    case DataType::Dense:
        return bedFile + ".ppbed";

    case DataType::SparseEigen:
        return bedFile +  ".eigen.sparsebed";

    case DataType::SparseRagged:
        return bedFile +  ".ragged.sparsebed";

    default:
        std::cerr << "ppFileForType - unsupported DataType: "
             << type
             << std::endl;
        assert(false);
        return {};
    }
}

std::string ppIndexFileForType(DataType type, const std::string &bedFile)
{
    switch (type) {
    case DataType::Dense:
        return bedFile +  ".ppbedindex";

    case DataType::SparseEigen:
        return bedFile +  ".eigen.sparsebedindex";

    case DataType::SparseRagged:
        return bedFile + ".ragged.sparsebedindex";

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
