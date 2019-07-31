#ifndef EIGENBAYESRKERNEL_H
#define EIGENBAYESRKERNEL_H

#include "sparsebayesrkernel.h"
#include "eigensparsemarker.h"

struct EigenBayesRKernel : public SparseBayesRKernel
{
    explicit EigenBayesRKernel(const EigenSparseMarker *marker);

    const VectorXd *ones = nullptr;

    VectorXdPtr calculateEpsilonChange(const double beta_old,
                                       const double beta) override;

protected:
    const EigenSparseMarker *esm = nullptr;

    double dot(const VectorXd &epsilon) const override;
};

#endif // EIGENBAYESRKERNEL_H
