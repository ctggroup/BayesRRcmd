#include "sparsebayesrkernel.h"

SparseBayesRKernel::SparseBayesRKernel(const std::shared_ptr<const SparseMarker> &marker)
    : BayesRKernel(marker)
    , sm(marker.get())
{
    assert(sm);
}

double SparseBayesRKernel::computeNum(const VectorXd &epsilon, const double beta_old)
{
    return computeNum(epsilon, beta_old, epsilonSum);
}

VectorXdPtr SparseBayesRKernel::calculateEpsilonChange(const double beta_old, const double beta)
{
    // now every update only saves delta epsilon sum
    epsilonSum = computeEpsilonSumUpdate(beta_old, beta);

    // return nullptr because we don't calculate the epsilon change here
    return nullptr;
}

double SparseBayesRKernel::computeNum(const VectorXd &epsilon, const double beta_old, const double epsilonSum)
{
    return beta_old * (static_cast<double>(sm->numInds) - 1.0) - sm->mean * epsilonSum / sm->sd + dot(epsilon);
}

double SparseBayesRKernel::computeEpsilonSumUpdate(const double beta_old, const double beta) const
{
    //Regardless of which scheme, the update of epsilonSum is the same
    const double dBeta = beta_old - beta;
    return dBeta * sm->Zsum / sm->sd - dBeta * sm->mean * static_cast<double>(sm->numInds) / sm->sd;
}
