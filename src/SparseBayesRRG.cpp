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

    sparseMarker->epsilonSum = m_isAsync ? m_asyncEpsilonSum : m_epsilonSum;
}

void SparseBayesRRG::writeWithUniqueLock(Marker *marker)
{
    auto* sparseMarker = dynamic_cast<SparseMarker*>(marker);
    assert(sparseMarker);

    if (m_isAsync)
        m_asyncEpsilonSum = sparseMarker->epsilonSum;
    else
        m_epsilonSum = sparseMarker->epsilonSum;
}

void SparseBayesRRG::updateGlobal(Marker *marker, const double beta_old, const double beta)
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    auto* sparseMarker = dynamic_cast<SparseMarker*>(marker);
    assert(sparseMarker);

    sparseMarker->updateEpsilon(m_epsilon, beta_old, beta);
}
