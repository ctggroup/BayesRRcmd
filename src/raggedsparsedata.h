#ifndef RAGGEDSPARSEDATA_H
#define RAGGEDSPARSEDATA_H

#include "sparsedata.h"

class RaggedSparseData : public SparseData
{
public:
    RaggedSparseData();

    using IndexVector = std::vector<int>;
    using RaggedVector = std::vector<IndexVector>;
    // vector containing the vectors the indexes of elements of the bed matrix which are one for each column
    RaggedVector Zones;
    // vector containing the vectors the indexes of elements of the bed matrix which are two for each column
    RaggedVector Ztwos;

    // vector containing the vectors the indexes of elements of the bed matrix which are missing for each column
    RaggedVector Zmissing;

    double dot(const unsigned int marker, const VectorXd &epsilon) const override;
    void updateEpsilon(VectorXd &epsilon, const unsigned int marker, const double beta_old, const double beta) const override;

protected:
    IndexVector* m_currentOnes = nullptr;
    IndexVector* m_currentTwos = nullptr;
    IndexVector* m_currentMissing = nullptr;

    void initialise() override;
    void beginSnpColumn(unsigned int snp) override;
    void processAllele(unsigned int snp, unsigned int individual, unsigned int allele1, unsigned int allele2) override;
    void endSnpColumn(unsigned int snp, unsigned int missingGenotypeCount) override;
};

#endif // RAGGEDSPARSEDATA_H
