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
