#ifndef DENSEBAYESWKERNEL_H
#define DENSEBAYESWKERNEL_H

#include "bayeswkernel.h"
#include "densemarker.h"

struct DenseBayesWKernel : public BayesWKernel
{
    explicit DenseBayesWKernel(const DenseMarker *marker);

    void setVi(const std::shared_ptr<VectorXd> &vi) override;
    void calculateSumFailure(const VectorXd &failure_vector);

    VectorXdPtr calculateEpsilonChange(const double beta) override;

    double exponent_sum() const override;
    double integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const override;

protected:
    const DenseMarker *dm = nullptr;
    std::shared_ptr<VectorXd> m_vi = nullptr;
};

#endif // DENSEBAYESWKERNEL_H
