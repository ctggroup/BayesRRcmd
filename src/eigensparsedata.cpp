#include "eigensparsedata.h"

const EigenSparseData::UnitDataType EigenSparseData::kMissing = std::numeric_limits<EigenSparseData::UnitDataType>::lowest();

EigenSparseData::EigenSparseData()
    : SparseData()
{

}

double EigenSparseData::dot(const unsigned int marker, const VectorXd &epsilon) const
{
    return Zg[marker].dot(epsilon) / sds(marker);
}

void EigenSparseData::updateEpsilon(VectorXd &epsilon, const unsigned int marker, const double beta_old, const double beta) const
{
    const double dBeta = beta_old - beta;
    epsilon += dBeta * Zg[marker] / sds(marker) - dBeta * means(marker) / sds(marker) * m_ones;
}

void EigenSparseData::initialise()
{
    Zg.resize(numSnps);

    // Data is expected to be about 80% zeros, so reserve a bit more than 20% of the expected space
    m_estimatedDataCount = static_cast<TupleList::size_type>(static_cast<double>(numInds) * 0.25);

    m_ones.setOnes(numInds);
}

void EigenSparseData::beginSnpColumn(unsigned int snp)
{
    (void) snp; // Unused

    // Create the triplet list for our sparse data representation
    m_tuples = std::make_shared<TupleList>();
    m_tuples->reserve(m_estimatedDataCount);
}

void EigenSparseData::processAllele(unsigned int snp, unsigned int individual, unsigned int allele1, unsigned int allele2)
{
    (void) snp; // Unused

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        m_tuples->emplace_back(individual, kMissing);
    } else if (allele1 == 1 || allele2 == 1) { // Not zero
        // Populate data for 1 or 2
        const double value = allele1 + allele2;
        m_tuples->emplace_back(individual, static_cast<UnitDataType>(value));
    }
}

void EigenSparseData::endSnpColumn(unsigned int snp, unsigned int missingGenotypeCount)
{
    const double mean = means[snp];

    // Update sqrdZ with imputed values
    sqrdZ[snp] += (mean * mean) * missingGenotypeCount;

    // Update Zsum with imputed values
    Zsum[snp] += mean * missingGenotypeCount;

    // Create the SparseVector
    auto &vector = Zg.at(snp);
    vector.resize(numInds); // Number of rows
    vector.reserve(m_tuples->size()); // Number of rows that are not zero

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
}
