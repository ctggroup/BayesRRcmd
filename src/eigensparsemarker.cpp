#include "eigensparsemarker.h"

void EigenSparseMarker::updateEpsilon(VectorXd &epsilon, const double beta_old, const double beta)
{
    SparseMarker::updateEpsilon(epsilon, beta_old, beta);

    const double dBeta = beta_old - beta;
    epsilon += dBeta * Zg / sd - dBeta * mean / sd * *ones;
}

double EigenSparseMarker::dot(const VectorXd &epsilon) const
{
    return Zg.dot(epsilon) / sd;
}
