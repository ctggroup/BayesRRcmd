#include "eigensparsemarkerbuilder.h"

static const EigenSparseMarker::UnitDataType kMissing = std::numeric_limits<EigenSparseMarker::UnitDataType>::lowest();

void EigenSparseMarkerBuilder::initialise(const unsigned int snp,
                                          const unsigned int numInds)
{
    MarkerBuilder::initialise(snp, numInds);

    m_marker.reset(new EigenSparseMarker);
    initialiseMarker();

    // Data is expected to be about 80% zeros, so reserve a bit more than 20% of the expected space
    const TupleList::size_type estimatedDataCount = static_cast<TupleList::size_type>(static_cast<double>(numInds) * 0.25);
    m_tuples = std::make_shared<TupleList>();
    m_tuples->reserve(estimatedDataCount);

    m_missingGenotypeCount = 0;
}

void EigenSparseMarkerBuilder::processAllele(unsigned int individual,
                                             unsigned int allele1,
                                             unsigned int allele2)
{
    auto* eigenMarker = dynamic_cast<EigenSparseMarker*>(m_marker.get());
    assert(eigenMarker);

    eigenMarker->updateStatistics(allele1, allele2);

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        m_tuples->emplace_back(individual, kMissing);
        ++m_missingGenotypeCount;
    } else if (allele1 == 1 || allele2 == 1) { // Not zero
        // Populate data for 1 or 2
        const double value = allele1 + allele2;
        m_tuples->emplace_back(individual, static_cast<UnitDataType>(value));
    }
}

void EigenSparseMarkerBuilder::endColumn()
{
    auto* eigenMarker = dynamic_cast<EigenSparseMarker*>(m_marker.get());
    assert(eigenMarker);

    // Calculate mean
    eigenMarker->mean /= m_numInds;
    const double mean = eigenMarker->mean;

    // Update sqrdZ with imputed values
    eigenMarker->sqrdZ += (mean * mean) * m_missingGenotypeCount;

    // Update Zsum with imputed values
    eigenMarker->Zsum += mean * m_missingGenotypeCount;

    // Calculate sd with imputed values
    eigenMarker->sd = std::sqrt((eigenMarker->sqrdZ - 2.0 * mean * eigenMarker->Zsum + m_numInds * mean * mean) /
                                 (m_numInds - 1.0));

    // Create the SparseVector
    auto &vector = eigenMarker->Zg;
    vector.resize(eigenMarker->numInds); // Number of rows
    vector.reserve(static_cast<Eigen::Index>(m_tuples->size())); // Number of rows that are not zero

    // Fill the SparseVector and doubles
    std::for_each(m_tuples->cbegin(), m_tuples->cend(), [&](const T &t) {
        const auto index = std::get<0>(t);
        auto value = std::get<1>(t);
        // If the value is zero or smaller, impute using the mean.
        if (value <= 0)
            value = mean;

        vector.insertBack(index) = value;
    });

    // Clean up
    m_tuples.reset();
    m_missingGenotypeCount = 0;
}
