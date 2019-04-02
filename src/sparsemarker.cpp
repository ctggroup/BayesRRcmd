#include "sparsemarker.h"

double SparseMarker::computeNum(VectorXd &epsilon, const double beta_old)
{
    return computeNum(epsilon, beta_old, epsilonSum);
}

void SparseMarker::updateEpsilon(VectorXd &epsilon, const double beta_old, const double beta)
{
    (void) epsilon; // Used by derived types
    epsilonSum += computeEpsilonSumUpdate(beta_old, beta);
}

double SparseMarker::computeNum(VectorXd &epsilon, const double beta_old, const double epsilonSum)
{
    return beta_old * (numInds - 1.0) - mean * epsilonSum / sd + dot(epsilon);
}

double SparseMarker::computeEpsilonSumUpdate(const double beta_old, const double beta) const
{
    //Regardless of which scheme, the update of epsilonSum is the same
    const double dBeta = beta_old - beta;
    return dBeta * Zsum / sd - dBeta * mean * numInds / sd;
}
