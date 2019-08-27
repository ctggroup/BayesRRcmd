#include "densebayesrkernel.h"

DenseRKernel::DenseRKernel(const std::shared_ptr<const DenseMarker> &marker)
    : BayesRKernel(marker)
    , dm(marker.get())
{
    assert(dm);
}

double DenseRKernel::computeNum(const VectorXd &epsilon, const double beta_old)
{
    //in order to not break async and sync updates for dense we change this
    //we now CX dot CX = N-1 given that both are already centered and scaled
    return dm->Cx->dot(epsilon) + beta_old * static_cast<double>(dm->numInds-1);
}

VectorXdPtr DenseRKernel::calculateEpsilonChange(const double beta_old, const double beta)
{
    return std::make_unique<VectorXd>((beta_old-beta) * *dm->Cx);
}
