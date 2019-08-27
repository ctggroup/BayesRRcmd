#ifndef RAGGEDBAYESRKERNEL_H
#define RAGGEDBAYESRKERNEL_H

#include "sparsebayesrkernel.h"
#include "raggedsparsemarker.h"

struct RaggedBayesRKernel : public SparseBayesRKernel
{
    explicit RaggedBayesRKernel(const std::shared_ptr<const RaggedSparseMarker> &marker);

    VectorXdPtr calculateEpsilonChange(const double beta_old,
                                       const double beta) override;

protected:
    const RaggedSparseMarker *rsm = nullptr;

    double dot(const VectorXd &epsilon) const override;
};

#endif // RAGGEDBAYESRKERNEL_H
