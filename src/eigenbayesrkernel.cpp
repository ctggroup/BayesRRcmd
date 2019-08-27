#include "eigenbayesrkernel.h"

EigenBayesRKernel::EigenBayesRKernel(const std::shared_ptr<const EigenSparseMarker> &marker)
    : SparseBayesRKernel(marker)
    , esm(marker.get())
{
    assert(esm);
}

VectorXdPtr EigenBayesRKernel::calculateEpsilonChange(const double beta_old, const double beta)
{
    SparseBayesRKernel::calculateEpsilonChange(beta_old, beta);
    const double dBeta = beta_old - beta;
    return std::make_unique<VectorXd>(dBeta * esm->Zg / esm->sd - dBeta * esm->mean / esm->sd * *ones);
}

double EigenBayesRKernel::dot(const VectorXd &epsilon) const
{
    return esm->Zg.dot(epsilon) / esm->sd;
}
