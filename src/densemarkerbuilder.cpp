#include "densemarkerbuilder.h"

#include "common.h"
#include "densemarker.h"

void DenseMarkerBuilder::initialise(const unsigned int snp,
                                    const unsigned int numInds)
{
    MarkerBuilder::initialise(snp, numInds);

    m_marker.reset(new DenseMarker);
    initialiseMarker();

    auto* denseMarker = dynamic_cast<DenseMarker*>(m_marker.get());
    assert(denseMarker);

    const unsigned int colSize = numInds * sizeof (double);
    denseMarker->buffer.reset(new unsigned char[colSize]);
    denseMarker->Cx = std::make_shared<Map<VectorXd>>(reinterpret_cast<double *>(denseMarker->buffer.get()),
                                                      numInds);
}

void DenseMarkerBuilder::processAllele(unsigned int individual,
                                       unsigned int allele1,
                                       unsigned int allele2)
{
    auto* denseMarker = dynamic_cast<DenseMarker*>(m_marker.get());
    assert(denseMarker);

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        // Don't store a marker value like this as it requires floating point comparisons later
        // which are not done properly. Instead, store the index of the individual in a vector and simply
        // iterate over the collected indices. Also means iterating over far fewer elements which may
        // make a noticeable difference as this scales up.
        m_missingIndices.push_back(individual);
    } else {
        const auto value = allele1 + allele2;
        (*denseMarker->Cx)[individual] = value;
        m_sum += value;
    }
}

void DenseMarkerBuilder::endColumn()
{
    auto* denseMarker = dynamic_cast<DenseMarker*>(m_marker.get());
    assert(denseMarker);

    auto& snpData = *denseMarker->Cx;

    const double mean = m_sum / (m_numInds - static_cast<double>(m_missingIndices.size()));
    std::for_each(m_missingIndices.cbegin(), m_missingIndices.cend(), [&](const unsigned int index) {
        snpData[index] = mean;
    });

    // Standardize genotypes
    snpData.array() -= snpData.mean();
    const auto sqn = snpData.squaredNorm();
    const auto sigma = 1.0 / (sqrt(sqn / (m_numInds - 1.0)));
    snpData.array() *= sigma;
}
