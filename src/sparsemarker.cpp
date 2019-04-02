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


void RaggedSparseMarker::updateEpsilon(VectorXd &epsilon, const double beta_old, const double beta)
{
    SparseMarker::updateEpsilon(epsilon, beta_old, beta);

    const double dBeta = beta_old - beta;
    const auto meanAdjustment = dBeta * mean / sd;
    // 1. Adjust for the means. If snp is 0, this will be the only adjustment made
    epsilon.array() -= meanAdjustment;

    // 2. Adjust for snp 1 values
    const double oneAdjustment = dBeta / sd;
    epsilon(Zones).array() += oneAdjustment;

    // 3. Adjust for snp 2 values
    epsilon(Ztwos).array() += 2 * oneAdjustment;

    // 4. For missing values, undo step 1
    epsilon(Zmissing).array() += meanAdjustment;
}

double RaggedSparseMarker::dot(const VectorXd &epsilon) const
{
    return (epsilon(Zones).sum() + 2 * epsilon(Ztwos).sum()) / sd;
}
