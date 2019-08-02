#include "raggedbayeswkernel.h"

RaggedBayesWKernel::RaggedBayesWKernel(const RaggedSparseMarker *marker)
    : BayesWKernel (marker)
    , rsm(marker)
{
    assert(marker);
}

void RaggedBayesWKernel::setVi(const VectorXd &vi)
{
    vi_sum = vi.sum();
    vi_2 = vi(rsm->Ztwos).sum();
    vi_1 = vi(rsm->Zones).sum();
    vi_0 = vi_sum - vi_1 - vi_2;
}

void RaggedBayesWKernel::calculateSumFailure(const VectorXd &failure_vector)
{
    int temp_sum = 0;
    temp_sum += failure_vector(rsm->Zones).sum();
    temp_sum += 2 * failure_vector(rsm->Ztwos).sum();

    sum_failure = (temp_sum - rsm->mean * failure_vector.array().sum()) / rsm->sd;
}

VectorXdPtr RaggedBayesWKernel::calculateEpsilonChange(const double beta)
{
    const auto mean_sd_ratio = rsm->mean / rsm->sd;
    const double meanAdjustment = mean_sd_ratio * beta;
    //Adjust for every memeber
    auto delta = std::make_unique<VectorXd>(VectorXd::Constant(rsm->numInds, -meanAdjustment));

    //And adjust even further for specific 1 and 2 allele values
    const double oneAdjustment = beta / rsm->sd;
    (*delta)(rsm->Zones).array() += oneAdjustment;

    const double twoAdjustment = 2 * oneAdjustment;
    (*delta)(rsm->Ztwos).array() += twoAdjustment;

    return delta;
}

double RaggedBayesWKernel::exponent_sum() const
{
    return (vi_1 * (1 - 2 * rsm->mean) + 4 * (1-rsm->mean) * vi_2 + vi_sum * rsm->mean * rsm->mean) /(rsm->sd*rsm->sd);
}

double RaggedBayesWKernel::integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const
{
    const auto mean_sd_ratio = rsm->mean / rsm->sd;
    const double temp = -alpha *s*sum_failure*sqrt_2Ck_sigmab +
            vi_sum - exp(alpha*mean_sd_ratio*s*sqrt_2Ck_sigmab) *
            (vi_0 + vi_1 * exp(-alpha * s*sqrt_2Ck_sigmab/rsm->sd) + vi_2* exp(-2 * alpha * s*sqrt_2Ck_sigmab/rsm->sd))
            -pow(s,2);
    return exp(temp);
}
