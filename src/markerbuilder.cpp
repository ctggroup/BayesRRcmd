#include "markerbuilder.h"

#include <iostream>

MarkerBuilder::~MarkerBuilder()
{

}

void MarkerBuilder::initialise(const unsigned int snp,
                               const unsigned int numInds)
{
    reset(numInds);
    m_snp = snp;
}

void MarkerBuilder::decompress(unsigned char *data, const IndexEntry &index) const
{
    assert(m_marker);
    m_marker->decompress(data, index);
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
    m_snp = 0;
    m_numInds = numInds;
    m_missingIndices.clear();
    m_sum = 0;
}

void MarkerBuilder::initialiseMarker()
{
    assert(m_marker);

    m_marker->i = m_snp;
    m_marker->numInds = static_cast<unsigned int>(m_numInds);
}
