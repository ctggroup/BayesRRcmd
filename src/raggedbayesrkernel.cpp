#include "raggedbayesrkernel.h"

RaggedBayesRKernel::RaggedBayesRKernel(const std::shared_ptr<const RaggedSparseMarker> &marker)
    : SparseBayesRKernel (marker)
    , rsm(marker.get())
{
    assert(rsm);
}

VectorXdPtr RaggedBayesRKernel::calculateEpsilonChange(const double beta_old, const double beta)
{
    SparseBayesRKernel::calculateEpsilonChange(beta_old, beta);

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

double RaggedBayesRKernel::dot(const VectorXd &epsilon) const
{
    return (epsilon(rsm->Zones).sum() + 2 * epsilon(rsm->Ztwos).sum()) / rsm->sd;
}
