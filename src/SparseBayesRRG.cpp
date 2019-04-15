#include "SparseBayesRRG.hpp"

#include "common.h"
#include "eigensparsemarker.h"
#include "sparsemarker.h"

SparseBayesRRG::SparseBayesRRG(const Data *data, Options &opt)
    : BayesRBase(data, opt)
{
}

SparseBayesRRG::~SparseBayesRRG()
{

}

MarkerBuilder *SparseBayesRRG::markerBuilder() const
{
    switch (m_opt.dataType) {
    case DataType::SparseEigen:
        // Fall through
    case DataType::SparseRagged:
        return builderForType(m_opt.dataType);

    default:
        std::cerr << "SparseBayesRRG::markerBuilder - unsupported type: "
                  << m_opt.dataType
                  << std::endl;
    }

    return nullptr;
}

void SparseBayesRRG::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    BayesRBase::init(K, markerCount, individualCount);

    m_asyncEpsilon = VectorXd(individualCount);

    m_asyncEpsilonSum = m_epsilonSum;

    m_ones.setOnes(individualCount);
}

void SparseBayesRRG::prepareForAnylsis()
{
    if (m_isAsync) {
        std::memcpy(m_asyncEpsilon.data(), m_epsilon.data(), static_cast<size_t>(m_epsilon.size()) * sizeof(double));
        m_asyncEpsilonSum = m_epsilonSum;
    }
}

void SparseBayesRRG::prepare(Marker *marker)
{
    // Hmmm
    if (auto* eigenSparseMarker = dynamic_cast<EigenSparseMarker*>(marker)) {
        eigenSparseMarker->ones = &m_ones;
    }
}

void SparseBayesRRG::readWithSharedLock(Marker *marker)
{
    auto* sparseMarker = dynamic_cast<SparseMarker*>(marker);
    assert(sparseMarker);
    //now we update to the global epsilonSum 
    sparseMarker->epsilonSum=m_epsilonSum;
}

void SparseBayesRRG::writeWithUniqueLock(Marker *marker)
{
    auto* sparseMarker = dynamic_cast<SparseMarker*>(marker);
    assert(sparseMarker);
    if (m_isAsync)
      {} //now the global node is in charge of updating m_epsilon  
    else
        m_epsilonSum += sparseMarker->epsilonSum;
}

void SparseBayesRRG::updateGlobal(Marker *marker, const double beta_old, const double beta,VectorXd& deltaEps)
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    auto* sparseMarker = dynamic_cast<SparseMarker*>(marker);
    assert(sparseMarker);
    m_epsilon+= deltaEps ;     //now epsilon=epsilon + 0 + update of epsilon. If vectorised this operation should not be expensive
    m_epsilonSum+=sparseMarker->epsilonSum; // now epsilonSum contains only deltaEpsilonSum
    m_betasqn+=beta*beta-beta_old*beta_old; // we move the squared norm computation to the global node, to avoid locking
}
