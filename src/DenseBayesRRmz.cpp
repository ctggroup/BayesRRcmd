/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "DenseBayesRRmz.hpp"
#include "densebayesrkernel.h"
#include "densemarker.h"
#include "common.h"

DenseBayesRRmz::DenseBayesRRmz(const Data *data, const Options *opt)
    : BayesRBase(data, opt)
{

}

DenseBayesRRmz::~DenseBayesRRmz()
{
}

std::unique_ptr<Kernel> DenseBayesRRmz::kernelForMarker(const ConstMarkerPtr &marker) const
{
    const auto denseMarker = dynamic_pointer_cast<const DenseMarker>(marker);
    assert(denseMarker);
    return std::make_unique<DenseRKernel>(denseMarker);
}

MarkerBuilder *DenseBayesRRmz::markerBuilder() const
{
    return builderForType(PreprocessDataType::Dense);
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

//update for mu, in dense case the cache variable m_epsilonSum is not used
void DenseBayesRRmz::updateMu(double old_mu,double N)
{
    m_epsilon = m_epsilon.array() + m_mu;// for dense and sparse we substract previous value
    m_mu = m_dist.norm_rng(m_epsilon.sum() / N, m_sigmaE / N); //update mu with the sum reduction 
    m_epsilon = m_epsilon.array() - m_mu;// for dense and sparse we substract again now epsilon =Y-mu-X*beta
   //we perform the equivalent update in epsilonSum for sparse this is important, for dense its ineffec.
}
