#include "SparseBayesRRG.hpp"

#include "common.h"
#include "eigenbayesrkernel.h"
#include "raggedbayesrkernel.h"
#include "sparsemarker.h"

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

SparseBayesRRG::SparseBayesRRG(const Data *data, const Options *opt)
    : BayesRBase(data, opt)
{
}

SparseBayesRRG::~SparseBayesRRG()
{

}

std::unique_ptr<Kernel> SparseBayesRRG::kernelForMarker(const ConstMarkerPtr &marker) const
{
    switch (m_opt->preprocessDataType) {
    case PreprocessDataType::SparseEigen:
    {
        const auto eigenSparseMarker = dynamic_pointer_cast<const EigenSparseMarker>(marker);
        assert(eigenSparseMarker);
        return std::make_unique<EigenBayesRKernel>(eigenSparseMarker);
    }

    case PreprocessDataType::SparseRagged:
    {
        const auto raggedSparseMarker = dynamic_pointer_cast<const RaggedSparseMarker>(marker);
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

    m_ones.setOnes(individualCount);
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

void SparseBayesRRG::resetAccumulators()
{
    BayesRBase::resetAccumulators();

    m_accumulatedEpsilonSum = 0;
}

void SparseBayesRRG::updateGlobal(const KernelPtr& kernel, const ConstAsyncResultPtr &result)
{
    BayesRBase::updateGlobal(kernel, result);

    auto* sparseKernel = dynamic_cast<SparseBayesRKernel*>(kernel.get());
    assert(sparseKernel);

    std::unique_lock lock(m_mutex);
    m_epsilonSum += sparseKernel->epsilonSum; // now epsilonSum contains only deltaEpsilonSum
}

void SparseBayesRRG::accumulate(const KernelPtr &kernel, const ConstAsyncResultPtr &result)
{
#ifdef MPI_TIMING_ENABLED
    const auto start = std::chrono::steady_clock::now();
#endif

    BayesRBase::accumulate(kernel, result);

    auto* sparseKernel = dynamic_cast<SparseBayesRKernel*>(kernel.get());
    assert(sparseKernel);

    std::unique_lock lock(m_accumulatorMutex);
    m_accumulatedEpsilonSum += sparseKernel->epsilonSum;

#ifdef MPI_TIMING_ENABLED
    const auto end = std::chrono::steady_clock::now();
    m_accumulateTime += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000000.0;
#endif
}

void SparseBayesRRG::updateMpi()
{
#ifdef MPI_ENABLED
#ifdef MPI_TIMING_ENABLED
    const auto start = std::chrono::steady_clock::now();
#endif
    // Take a copy of the accumulated values
    const auto localEpsilonSum = m_accumulatedEpsilonSum;

    // MPI_Allreduce
    MPI_Allreduce(&localEpsilonSum, &m_accumulatedEpsilonSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifdef MPI_TIMING_ENABLED
    const auto mpiSync = std::chrono::steady_clock::now();
#endif

    // Subtract local accumulations for global
    m_accumulatedEpsilonSum -= localEpsilonSum;

    // Apply accumulations from other processes
    m_epsilonSum += m_accumulatedEpsilonSum;

    // Call last - it calls resetAccumulators
    BayesRBase::updateMpi();

#ifdef MPI_TIMING_ENABLED
    const auto end = std::chrono::steady_clock::now();

    const auto mpiCount = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(mpiSync - start).count()) / 1000000.0;

    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    std::vector<double> waitTime(worldSize);
    MPI_Gather(&mpiCount, 1, MPI_DOUBLE, waitTime.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::transform(m_waitTime.begin(), m_waitTime.end(), waitTime.begin(),
                   m_waitTime.begin(), std::plus<double>());

    const auto totalCount = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000000.0;

    std::vector<double> mpiTime(worldSize);
    MPI_Gather(&totalCount, 1, MPI_DOUBLE, mpiTime.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::transform(m_mpiTime.begin(), m_mpiTime.end(), mpiTime.begin(),
                   m_mpiTime.begin(), std::plus<double>());
#endif
#endif
}


void SparseBayesRRG::updateMu(double old_mu,double N)
{
    m_epsilon = m_epsilon.array() + m_mu;// for dense and sparse we substract previous value
    m_epsilonSum+=old_mu*double(N); //for sparse this is important for dense its ineffectual
    m_mu = m_dist.norm_rng(m_epsilonSum / N, m_sigmaE / N); //update mu with the sum reduction 
    m_epsilon = m_epsilon.array() - m_mu;// for dense and sparse we substract again now epsilon =Y-mu-X*beta
    m_epsilonSum-=m_mu*N;//we perform the equivalent update in epsilonSum for sparse this is important, for dense its ineffec.
}
