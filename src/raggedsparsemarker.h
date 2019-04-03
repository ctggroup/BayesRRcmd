#ifndef RAGGEDSPARSEMARKER_H
#define RAGGEDSPARSEMARKER_H

#include "sparsemarker.h"
#include "markerbuilder.h"

struct RaggedSparseMarker : public SparseMarker
{
    using IndexVector = std::vector<int>;
    // the indexes of elements of the bed matrix which are one for this column
    IndexVector Zones;
    // the indexes of elements of the bed matrix which are two for this column
    IndexVector Ztwos;

    // the indexes of elements of the bed matrix which are missing for this column
    IndexVector Zmissing;

    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;

protected:
    double dot(const VectorXd &epsilon) const override;
};


template<>
CompressedMarker compress(const RaggedSparseMarker* marker);


template<>
void write(const RaggedSparseMarker* marker, std::ostream *outStream);

template<>
void MarkerBuilder<RaggedSparseMarker>::initialise(const unsigned int snp,
                                                   const unsigned int numInds);

template<>
void MarkerBuilder<RaggedSparseMarker>::processAllele(unsigned int individual,
                                                      unsigned int allele1,
                                                      unsigned int allele2);
template<>
void MarkerBuilder<RaggedSparseMarker>::endColumn();

#endif // RAGGEDSPARSEMARKER_H
