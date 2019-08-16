#ifndef BAYESWKERNEL_H
#define BAYESWKERNEL_H

#include "kernel.h"

#include <Eigen/Eigen>

using namespace Eigen;

struct BayesWKernel : public Kernel
{
    explicit BayesWKernel(const ConstMarkerPtr &marker) : Kernel(marker) {}
    ~BayesWKernel();


    double sum_failure = 0;

    virtual void setVi(const std::shared_ptr<VectorXd>& vi) = 0;
    // Should really be done as part of the preprocess step
    virtual void calculateSumFailure(const VectorXd &failure_vector) = 0;

    virtual VectorXdPtr calculateResidualUpdate(const double beta) = 0;
    virtual VectorXdPtr calculateEpsilonChange(const double beta_old, const double beta) = 0;

    virtual double exponent_sum() const = 0;
    virtual double integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const = 0;

protected:
};

#endif // BAYESWKERNEL_H
