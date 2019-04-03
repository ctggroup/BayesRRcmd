#include "markerbuilder.h"

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
    reset();
    return m_marker.release();
}

void MarkerBuilder::reset(const double numInds)
{
    m_numInds = numInds;
    m_missingIndices.clear();
    m_sum = 0;
}
