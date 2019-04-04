#include "common.h"

#include "densemarkerbuilder.h"
#include "raggedsparsemarkerbuilder.h"

#include <cassert>
#include <iostream>

MarkerBuilder *builderForType(const DataType type)
{
    switch (type)
    {
    case DataType::Dense:
        return new DenseMarkerBuilder;

    case DataType::SparseRagged:
        return new RaggedSparseMarkerBuilder;

    case DataType::None:
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

    case DataType::SparseRagged:
        return bedFile +  ".ragged.sparsebed";

    default:
        std::cerr << "ppFileForType - unsupported DataType: "
             << type
             << std::endl;
        return {};
    }
}

std::string ppIndexFileForType(DataType type, const std::string &bedFile)
{
    switch (type) {
    case DataType::Dense:
        return bedFile +  ".ppbedindex";

    case DataType::SparseRagged:
        return bedFile + ".ragged.sparsebedindex";

    default:
        std::cerr << "ppIndexFileForType - unsupported DataType: "
             << type
             << std::endl;
        return {};
    }
}
