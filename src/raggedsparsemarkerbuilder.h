#ifndef RAGGEDMARKERBUILDER_H
#define RAGGEDMARKERBUILDER_H

#include "markerbuilder.h"

class RaggedSparseMarkerBuilder : public MarkerBuilder
{
public:
    explicit RaggedSparseMarkerBuilder() = default;

    void initialise(const unsigned int snp,
                    const unsigned int numInds) override;

    void processAllele(unsigned int individual,
                       unsigned int allele1,
                       unsigned int allele2) override;

    void endColumn() override;
};

#endif // RAGGEDMARKERBUILDER_H
