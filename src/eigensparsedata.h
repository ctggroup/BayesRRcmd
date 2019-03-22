#ifndef EIGENSPARSEDATA_H
#define EIGENSPARSEDATA_H

#include "sparsedata.h"

class EigenSparseData : public SparseData
{
protected:
    using UnitDataType = double;
    using T = std::tuple<Eigen::Index, UnitDataType>;
    using TupleList = std::vector<T>;

public:
    EigenSparseData();

    // this is a sparse matrix that contains the uncentered and unscaled elements of the bed matrix
    std::vector<SparseVector<UnitDataType>> Zg;

    double dot(const unsigned int marker, const VectorXd &epsilon) const override;
    void updateEpsilon(VectorXd &epsilon, const unsigned int marker, const double beta_old, const double beta) const override;

protected:
    static const UnitDataType kMissing;

    TupleList::size_type m_estimatedDataCount = 0;

    using TupleListPtr = std::shared_ptr<TupleList>;
    TupleListPtr m_tuples = nullptr;

    VectorXd m_ones;

    void initialise() override;
    void beginSnpColumn(unsigned int snp) override;
    void processAllele(unsigned int snp, unsigned int individual, unsigned int allele1, unsigned int allele2) override;
    void endSnpColumn(unsigned int snp, unsigned int missingGenotypeCount) override;
};

#endif // EIGENSPARSEDATA_H
