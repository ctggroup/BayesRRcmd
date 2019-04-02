#ifndef MARKERBUILDER_H
#define MARKERBUILDER_H

#include <memory>
#include <vector>

template <typename MarkerType>
class MarkerBuilder
{
public:
    void initialise(const unsigned int snp, const double numInds);
    void processAllele(unsigned int individual,
                       unsigned int allele1,
                       unsigned int allele2);
    void endColumn();
    MarkerType* build();

protected:
    std::unique_ptr<MarkerType> m_marker = nullptr;

    double m_numInds = 0;
    std::vector<unsigned int> m_missingIndices;
    double m_sum = 0;

    void reset(const double numInds = 0);
};

template<typename MarkerType>
MarkerType *MarkerBuilder<MarkerType>::build()
{
    reset();
    return m_marker.release();
}

template<typename MarkerType>
void MarkerBuilder<MarkerType>::reset(const double numInds)
{
    m_numInds = numInds;
    m_missingIndices.clear();
    m_sum = 0;
}

#endif // MARKERBUILDER_H
