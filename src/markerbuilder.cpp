#include "markerbuilder.h"

#include <iostream>

MarkerBuilder::~MarkerBuilder()
{

}

void MarkerBuilder::initialise(const unsigned int snp,
                               const unsigned int numInds)
{
    (void) snp; // Unused
    reset(numInds);
}

Marker *MarkerBuilder::build()
{
    if (!m_marker) {
        std::cerr << "Not building Marker, was MarkerBuilder::initialise called?" << std::endl;
        return nullptr;
    }

    if (!m_marker->isValid()) {
        std::cerr << "Marker is not valid!" << std::endl;
        return nullptr;
    }

    reset();
    return m_marker.release();
}

void MarkerBuilder::reset(const double numInds)
{
    m_numInds = numInds;
    m_missingIndices.clear();
    m_sum = 0;
}
