#ifndef RAGGEDBAYESWKERNEL_H
#define RAGGEDBAYESWKERNEL_H

#include "bayeswkernel.h"
#include "raggedsparsemarker.h"

struct RaggedBayesWKernel : public BayesWKernel
{
    explicit RaggedBayesWKernel(const std::shared_ptr<const RaggedSparseMarker> &marker);

    double vi_sum = 0;
    double vi_2 = 0;
    double vi_1 = 0;
    double vi_0 = 0;

    void setVi(const std::shared_ptr<VectorXd> &vi) override;
    void calculateSumFailure(const VectorXd &failure_vector);

    VectorXdPtr calculateResidualUpdate(const double beta) override;
    VectorXdPtr calculateEpsilonChange(const double beta_old, const double beta) override;

    double exponent_sum() const override;
    double integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const override;

protected:
    const RaggedSparseMarker *rsm = nullptr;
};

#endif // RAGGEDBAYESWKERNEL_H
