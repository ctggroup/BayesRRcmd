#include "raggedsparsedata.h"

RaggedSparseData::RaggedSparseData()
    : SparseData()
{

}

double RaggedSparseData::dot(const unsigned int marker, const VectorXd &epsilon) const
{
    return (epsilon(Zones[marker]).sum() + 2 * epsilon(Ztwos[marker]).sum()) / sds(marker);
}

void RaggedSparseData::updateEpsilon(VectorXd &epsilon, const unsigned int marker, const double beta_old, const double beta) const
{
    const double dBeta = beta_old - beta;
    const auto meanAdjustment = dBeta * means(marker) / sds(marker);
    // 1. Adjust for the means. If snp is 0, this will be the only adjustment made
    epsilon.array() -= meanAdjustment;

    // 2. Adjust for snp 1 values
    const double oneAdjustment = dBeta / sds(marker);
    epsilon(Zones[marker]).array() += oneAdjustment;

    // 3. Adjust for snp 2 values
    epsilon(Ztwos[marker]).array() += 2 * oneAdjustment;

    // 4. For missing values, undo step 1
    epsilon(Zmissing[marker]).array() += meanAdjustment;
}

void RaggedSparseData::initialise()
{
    // Initialise our RaggedVectors with IndexVectors with the estimated amount of reserved space.
    using size_type = RaggedVector::size_type;
    size_type estimatedDataCount = static_cast<size_type>(static_cast<double>(numInds) * 0.2);
    IndexVector tmp;
    tmp.reserve(estimatedDataCount);

    Zones = RaggedVector(numSnps, tmp);
    Ztwos = RaggedVector(numSnps, tmp);

    size_type estimatedMissingCount = static_cast<size_type>(static_cast<double>(numInds) * 0.01);
    tmp.reserve(estimatedMissingCount);

    Zmissing = RaggedVector(numSnps, tmp);
}

void RaggedSparseData::beginSnpColumn(unsigned int snp)
{
    m_currentOnes = &Zones[snp];
    m_currentTwos = &Ztwos[snp];
    m_currentMissing = &Zmissing[snp];
}

void RaggedSparseData::processAllele(unsigned int snp, unsigned int individual, unsigned int allele1, unsigned int allele2)
{
    assert(m_currentOnes);
    assert(m_currentTwos);
    assert(m_currentMissing);

    (void) snp; // Unused

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        m_currentMissing->emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 0) {
        m_currentOnes->emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 1) {
        m_currentTwos->emplace_back(individual);
    }
}

void RaggedSparseData::endSnpColumn(unsigned int snp, unsigned int missingGenotypeCount)
{
    (void) snp; // Unused
    (void) missingGenotypeCount; // Unused

    // No need to impute missing values

    m_currentOnes = nullptr;
    m_currentTwos = nullptr;
    m_currentMissing = nullptr;
}
