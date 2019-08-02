#include "SparseBayesRRG.hpp"

#include "common.h"
#include "eigenbayesrkernel.h"
#include "raggedbayesrkernel.h"
#include "sparsemarker.h"

SparseBayesRRG::SparseBayesRRG(const Data *data, const Options *opt)
    : BayesRBase(data, opt)
{
}

SparseBayesRRG::~SparseBayesRRG()
{

}

std::unique_ptr<Kernel> SparseBayesRRG::kernelForMarker(const Marker *marker) const
{
    switch (m_opt->preprocessDataType) {
    case PreprocessDataType::SparseEigen:
    {
        const auto* eigenSparseMarker = dynamic_cast<const EigenSparseMarker*>(marker);
        assert(eigenSparseMarker);
        return std::make_unique<EigenBayesRKernel>(eigenSparseMarker);
    }

    case PreprocessDataType::SparseRagged:
    {
        const auto* raggedSparseMarker = dynamic_cast<const RaggedSparseMarker*>(marker);
        assert(raggedSparseMarker);
        return std::make_unique<RaggedBayesRKernel>(raggedSparseMarker);
    }

    default:
        std::cerr << "SparseBayesRRG::kernelForMarker - unsupported type: "
                  << m_opt->preprocessDataType
                  << std::endl;
    }

    return {};
}

MarkerBuilder *SparseBayesRRG::markerBuilder() const
{
    switch (m_opt->preprocessDataType) {
    case PreprocessDataType::SparseEigen:
        // Fall through
    case PreprocessDataType::SparseRagged:
        return builderForType(m_opt->preprocessDataType);

    default:
        std::cerr << "SparseBayesRRG::markerBuilder - unsupported type: "
                  << m_opt->preprocessDataType
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

void SparseBayesRRG::prepare(BayesRKernel *kernel)
{
    // Hmmm
    if (auto* eigenBayesRKernel = dynamic_cast<EigenBayesRKernel*>(kernel)) {
        eigenBayesRKernel->ones = &m_ones;
    }
}

void SparseBayesRRG::readWithSharedLock(BayesRKernel *kernel)
{
    auto* sparseKernel = dynamic_cast<SparseBayesRKernel*>(kernel);
    assert(sparseKernel);
    //now we update to the global epsilonSum 
    sparseKernel->epsilonSum = m_epsilonSum;
}

void SparseBayesRRG::writeWithUniqueLock(BayesRKernel *kernel)
{
    auto* sparseKernel = dynamic_cast<SparseBayesRKernel*>(kernel);
    assert(sparseKernel);
    if (m_isAsync)
      {} //now the global node is in charge of updating m_epsilon  
    else
        m_epsilonSum += sparseKernel->epsilonSum;
}

void SparseBayesRRG::updateGlobal(Kernel *kernel, const double beta_old, const double beta, const VectorXd& deltaEps)
{
    BayesRBase::updateGlobal(kernel, beta_old, beta, deltaEps);

    auto* sparseKernel = dynamic_cast<SparseBayesRKernel*>(kernel);
    assert(sparseKernel);
    m_epsilonSum += sparseKernel->epsilonSum; // now epsilonSum contains only deltaEpsilonSum
}


void SparseBayesRRG::updateMu(double old_mu,double N)
{
    m_epsilon = m_epsilon.array() + m_mu;// for dense and sparse we substract previous value
    m_epsilonSum+=old_mu*double(N); //for sparse this is important for dense its ineffectual
    m_mu = m_dist.norm_rng(m_epsilonSum / N, m_sigmaE / N); //update mu with the sum reduction 
    m_epsilon = m_epsilon.array() - m_mu;// for dense and sparse we substract again now epsilon =Y-mu-X*beta
    m_epsilonSum-=m_mu*N;//we perform the equivalent update in epsilonSum for sparse this is important, for dense its ineffec.
}
