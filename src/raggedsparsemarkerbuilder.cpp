#include "raggedsparsemarkerbuilder.h"

#include "raggedsparsemarker.h"

RaggedSparseMarkerBuilder::RaggedSparseMarkerBuilder()
    : MarkerBuilder ()
{

}

void RaggedSparseMarkerBuilder::initialise(const unsigned int snp,
                                           const unsigned int numInds)
{
    MarkerBuilder::initialise(snp, numInds);

    m_marker.reset(new RaggedSparseMarker);
    m_marker->i = snp;
    m_marker->numInds = numInds;

    auto* raggedMarker = dynamic_cast<RaggedSparseMarker*>(m_marker.get());
    assert(raggedMarker);

    using size_type = RaggedSparseMarker::IndexVector::size_type;
    const size_type estimatedDataCount = static_cast<size_type>(numInds * 0.2);
    raggedMarker->Zones.reserve(estimatedDataCount);
    raggedMarker->Ztwos.reserve(estimatedDataCount);

    const size_type estimatedMissingCount = static_cast<size_type>(numInds * 0.01);
    raggedMarker->Zmissing.reserve(estimatedMissingCount);
}

void RaggedSparseMarkerBuilder::processAllele(unsigned int individual,
                                              unsigned int allele1,
                                              unsigned int allele2)
{
    auto* raggedMarker = dynamic_cast<RaggedSparseMarker*>(m_marker.get());
    assert(raggedMarker);

    raggedMarker->updateStatistics(allele1, allele2);

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        raggedMarker->Zmissing.emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 0) {
        raggedMarker->Zones.emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 1) {
        raggedMarker->Ztwos.emplace_back(individual);
    }
}

void RaggedSparseMarkerBuilder::endColumn()
{
    auto* raggedMarker = dynamic_cast<RaggedSparseMarker*>(m_marker.get());
    assert(raggedMarker);

    // Calculate mean
    raggedMarker->mean /= static_cast<double>(raggedMarker->Zmissing.size());

    // Calculate sd
    const double mean = raggedMarker->mean;
    raggedMarker->sd = std::sqrt((raggedMarker->sqrdZ - 2.0 * mean * raggedMarker->Zsum + m_numInds * mean * mean) /
                                 (m_numInds - 1.0));
}

void RaggedSparseMarkerBuilder::decompress(unsigned char *data,
                                           const IndexEntry &index) const
{
    auto* raggedMarker = dynamic_cast<RaggedSparseMarker*>(m_marker.get());
    assert(raggedMarker);

    raggedMarker->decompress(data, index);
}
