#ifndef RAGGEDBAYESWKERNEL_H
#define RAGGEDBAYESWKERNEL_H

#include "bayeswkernel.h"
#include "raggedsparsemarker.h"

struct RaggedBayesWKernel : public BayesWKernel
{
    explicit RaggedBayesWKernel(const RaggedSparseMarker *marker);

    double vi_sum = 0;
    double vi_2 = 0;
    double vi_1 = 0;
    double vi_0 = 0;

    void setVi(const VectorXd &vi) override;
    void calculateSumFailure(const VectorXd &failure_vector);

    VectorXdPtr calculateEpsilonChange(const double beta) override;

    double exponent_sum() const override;
    double integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const override;

protected:
    const RaggedSparseMarker *rsm = nullptr;
};

#endif // RAGGEDBAYESWKERNEL_H
