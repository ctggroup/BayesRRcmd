#ifndef BAYESRKERNEL_H
#define BAYESRKERNEL_H

#include "kernel.h"

struct BayesRKernel : public Kernel
{
    explicit BayesRKernel(const ConstMarkerPtr &marker) : Kernel(marker) {}
    ~BayesRKernel();

    virtual double computeNum(const VectorXd &epsilon,
                              const double beta_old) = 0;

    virtual VectorXdPtr calculateEpsilonChange(const double beta_old,
                                               const double beta) = 0;
};

#endif // BAYESRKERNEL_H
