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

    CompressedMarker compress() const override;

    size_t size() const override;
    void write(std::ostream *outStream) const override;

protected:
    double dot(const VectorXd &epsilon) const override;
};

#endif // RAGGEDSPARSEMARKER_H
