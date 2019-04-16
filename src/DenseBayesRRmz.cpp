/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "DenseBayesRRmz.hpp"
#include "densemarker.h"
#include "common.h"

DenseBayesRRmz::DenseBayesRRmz(const Data *data, Options &opt)
    : BayesRBase (data, opt)
{

}

DenseBayesRRmz::~DenseBayesRRmz()
{
}

MarkerBuilder *DenseBayesRRmz::markerBuilder() const
{
    return builderForType(DataType::Dense);
}

void DenseBayesRRmz::updateGlobal(Marker *marker, const double beta_old, const double beta,VectorXd& deltaEps )
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    auto* denseMarker = dynamic_cast<DenseMarker*>(marker);
    assert(denseMarker);

    m_epsilon += deltaEps;
    m_betasqn+=beta*beta-beta_old*beta_old;
}

void DenseBayesRRmz::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    BayesRBase::init(K, markerCount, individualCount);

    if (m_isAsync)
        m_asyncEpsilon = VectorXd(individualCount);
}

void DenseBayesRRmz::prepareForAnylsis()
{
    if (m_isAsync)
        std::memcpy(m_asyncEpsilon.data(), m_epsilon.data(), static_cast<size_t>(m_epsilon.size()) * sizeof(double));
}

void DenseBayesRRmz::readWithSharedLock(Marker *marker)
{
    auto* denseMarker = dynamic_cast<DenseMarker*>(marker);
    assert(denseMarker);

    denseMarker->component = m_components(denseMarker->i);
}

//update for mu, in dense case the cache variable m_epsilonSum is not used
void DenseBayesRRmz::updateMu(double old_mu,double N)
{
    m_epsilon = m_epsilon.array() + m_mu;// for dense and sparse we substract previous value
    m_mu = m_dist.norm_rng(m_epsilon.sum() / N, m_sigmaE / N); //update mu with the sum reduction 
    m_epsilon = m_epsilon.array() - m_mu;// for dense and sparse we substract again now epsilon =Y-mu-X*beta
   //we perform the equivalent update in epsilonSum for sparse this is important, for dense its ineffec.
}
