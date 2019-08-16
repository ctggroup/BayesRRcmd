#ifndef DENSEBAYESRKERNEL_H
#define DENSEBAYESRKERNEL_H

#include "bayesrkernel.h"
#include "densemarker.h"

struct DenseRKernel : public BayesRKernel
{
    explicit DenseRKernel(const std::shared_ptr<const DenseMarker> &marker);

    double computeNum(const VectorXd &epsilon, const double beta_old) override;
    VectorXdPtr calculateEpsilonChange(const double beta_old, const double beta) override;

protected:
    const DenseMarker *dm = nullptr;
};

#endif // DENSEBAYESRKERNEL_H
