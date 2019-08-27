#ifndef EIGENSPARSEMARKERBUILDER_H
#define EIGENSPARSEMARKERBUILDER_H

#include "markerbuilder.h"

#include "eigensparsemarker.h"

class EigenSparseMarkerBuilder : public MarkerBuilder
{
public:
    explicit EigenSparseMarkerBuilder() = default;

    void initialise(const unsigned int snp,
                    const unsigned int numInds) override;

    void processAllele(unsigned int individual,
                       unsigned int allele1,
                       unsigned int allele2) override;

    void endColumn() override;

protected:
    using UnitDataType = EigenSparseMarker::UnitDataType;
    using T = std::tuple<Eigen::Index, UnitDataType>;
    using TupleList = std::vector<T>;
    using TupleListPtr = std::shared_ptr<TupleList>;
    TupleListPtr m_tuples = nullptr;

    double m_missingGenotypeCount = 0;
};

#endif // EIGENSPARSEMARKERBUILDER_H
