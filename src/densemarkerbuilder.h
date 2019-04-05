#ifndef DENSEMARKERBUILDER_H
#define DENSEMARKERBUILDER_H

#include "markerbuilder.h"

class DenseMarkerBuilder : public MarkerBuilder
{
public:
    explicit DenseMarkerBuilder() = default;

    void initialise(const unsigned int snp,
                    const unsigned int numInds) override;

    void processAllele(unsigned int individual,
                       unsigned int allele1,
                       unsigned int allele2) override;

    void endColumn() override;
};

#endif // DENSEMARKERBUILDER_H
