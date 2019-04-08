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

void DenseBayesRRmz::updateGlobal(Marker *marker, const double beta_old, const double beta)
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    auto* denseMarker = dynamic_cast<DenseMarker*>(marker);
    assert(denseMarker);

    m_epsilon -= *denseMarker->Cx * (beta - beta_old);
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
