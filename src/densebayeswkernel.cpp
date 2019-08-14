#include "densebayeswkernel.h"

DenseBayesWKernel::DenseBayesWKernel(const std::shared_ptr<const DenseMarker> &marker)
    : BayesWKernel(marker)
    , dm(marker.get())
{
    assert(dm);
}

void DenseBayesWKernel::setVi(const std::shared_ptr<VectorXd> &vi)
{
    m_vi = vi;
}

void DenseBayesWKernel::calculateSumFailure(const VectorXd &failure_vector)
{
    sum_failure = (dm->Cx->array() * failure_vector.array()).sum();
}

VectorXdPtr DenseBayesWKernel::calculateResidualUpdate(const double beta)
{
    return std::make_unique<VectorXd>(*dm->Cx * beta);
}

VectorXdPtr DenseBayesWKernel::calculateEpsilonChange(const double beta_old, const double beta)
{
    return std::make_unique<VectorXd>((beta_old-beta) * *dm->Cx);
}

double DenseBayesWKernel::exponent_sum() const
{
    assert(m_vi);
    return (m_vi->array() * dm->Cx->array() * dm->Cx->array()).sum();
}

double DenseBayesWKernel::integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const
{
    assert(m_vi);
    //vi is a vector of exp(vi)
    double temp = -alpha *s*sum_failure*sqrt_2Ck_sigmab + (m_vi->array()* (1 - (-dm->Cx->array()*s*sqrt_2Ck_sigmab*alpha).exp() )).sum() -pow(s,2);
    return exp(temp);
}
