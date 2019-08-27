#ifndef SPARSEBAYESRKERNEL_H
#define SPARSEBAYESRKERNEL_H

#include "bayesrkernel.h"
#include "sparsemarker.h"

struct SparseBayesRKernel : public BayesRKernel
{
    explicit SparseBayesRKernel(const std::shared_ptr<const SparseMarker> &marker);

    double epsilonSum = 0;

    double computeNum(const VectorXd &epsilon, const double beta_old) override;
    VectorXdPtr calculateEpsilonChange(const double beta_old, const double beta) override;

protected:
    const SparseMarker *sm = nullptr;

    virtual double computeNum(const VectorXd &epsilon,
                              const double beta_old,
                              const double epsilonSum);

    virtual double dot(const VectorXd &epsilon) const = 0;

    virtual double computeEpsilonSumUpdate(const double beta_old,
                                           const double beta) const;
};

#endif // SPARSEBAYESRKERNEL_H
