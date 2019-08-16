#include "raggedbayeswkernel.h"

RaggedBayesWKernel::RaggedBayesWKernel(const std::shared_ptr<const RaggedSparseMarker> &marker)
    : BayesWKernel (marker)
    , rsm(marker.get())
{
    assert(marker);
}

void RaggedBayesWKernel::setVi(const std::shared_ptr<VectorXd> &vi)
{
    vi_sum = vi->sum();
    vi_2 = (*vi)(rsm->Ztwos).sum();
    vi_1 = (*vi)(rsm->Zones).sum();
    vi_0 = vi_sum - vi_1 - vi_2;
}

void RaggedBayesWKernel::calculateSumFailure(const VectorXd &failure_vector)
{
    int temp_sum = 0;
    temp_sum += failure_vector(rsm->Zones).sum();
    temp_sum += 2 * failure_vector(rsm->Ztwos).sum();

    sum_failure = (temp_sum - rsm->mean * failure_vector.array().sum()) / rsm->sd;
}

VectorXdPtr RaggedBayesWKernel::calculateResidualUpdate(const double beta)
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

VectorXdPtr RaggedBayesWKernel::calculateEpsilonChange(const double beta_old, const double beta)
{
    const double dBeta = beta_old - beta;
    const auto meanAdjustment = dBeta * rsm->mean / rsm->sd;
    // 1. Adjust for the means. If snp is 0, this will be the only adjustment made
    auto delta = std::make_unique<VectorXd>(VectorXd::Constant(rsm->numInds, -meanAdjustment));

    // 2. Adjust for snp 1 values
    const double oneAdjustment = dBeta / rsm->sd;
    (*delta)(rsm->Zones).array() += oneAdjustment;

    // 3. Adjust for snp 2 values
    (*delta)(rsm->Ztwos).array() += 2 * oneAdjustment;

    // 4. For missing values, undo step 1
    (*delta)(rsm->Zmissing).array() += meanAdjustment;

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
