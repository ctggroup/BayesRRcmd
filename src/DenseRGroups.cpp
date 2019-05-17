
#include "DenseRGroups.hpp"
#include "densemarker.h"
#include "common.h"

DenseRGroups::DenseRGroups(const Data *data, Options &opt)
    : BayesRGroupsBase (data, opt)
{

}

DenseRGroups::~DenseRGroups()
{
}

MarkerBuilder *DenseRGroups::markerBuilder() const
{
   return builderForType(DataType::Dense);
}

void DenseRGroups::init(int K, unsigned int markerCount, unsigned int individualCount, unsigned int groupCount)
{
    BayesRGroupsBase::init(K, markerCount, individualCount, groupCount);

    if (m_isAsync)
        m_asyncEpsilon = VectorXd(individualCount);
}

void DenseRGroups::prepareForAnylsis()
{
    if (m_isAsync)
        std::memcpy(m_asyncEpsilon.data(), m_epsilon.data(), static_cast<size_t>(m_epsilon.size()) * sizeof(double));
}

void DenseRGroups::readWithSharedLock(Marker *marker)
{
    auto* denseMarker = dynamic_cast<DenseMarker*>(marker);
    assert(denseMarker);

    denseMarker->component = m_components(denseMarker->i);
}

//update for mu, in dense case the cache variable m_epsilonSum is not used
void DenseRGroups::updateMu(double old_mu,double N)
{
    m_epsilon = m_epsilon.array() + m_mu;// for dense and sparse we substract previous value
    m_mu = m_dist.norm_rng(m_epsilon.sum() / N, m_sigmaE / N); //update mu with the sum reduction
    m_epsilon = m_epsilon.array() - m_mu;// for dense and sparse we substract again now epsilon =Y-mu-X*beta
   //we perform the equivalent update in epsilonSum for sparse this is important, for dense its ineffec.
}


void DenseRGroups::updateGlobal(Marker *marker, const double beta_old, const double beta,VectorXd& deltaEps )
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    auto* denseMarker = dynamic_cast<DenseMarker*>(marker);
    assert(denseMarker);

    m_epsilon += deltaEps;
    m_betasqn+=beta*beta-beta_old*beta_old;
}
