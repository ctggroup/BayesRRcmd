#include "analysis.h"

#include <numeric>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

Analysis::Analysis(const Data *data, const Options *opt)
    : m_data(data)
    , m_opt(opt)
{
#ifdef MPI_ENABLED
    if (m_opt->useHybridMpi) {
        MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &m_worldSize);
    }
#endif
}

Analysis::~Analysis() = default;

IndexEntry Analysis::indexEntry(unsigned int i) const
{
    if (!m_data)
        return {};

    return m_data->ppbedIndex[i];
}

bool Analysis::compressed() const
{
    return m_opt->compress;
}

unsigned char* Analysis::compressedData() const
{
    if (!m_data)
        return nullptr;

    return reinterpret_cast<unsigned char*>(m_data->ppBedMap);
}

std::string Analysis::preprocessedFile() const
{
    return ppFileForType(*m_opt);
}

int Analysis::runGibbs(AnalysisGraph *analysis)
{
    std::vector<unsigned int> markerI(m_data->numSnps);
    std::iota(markerI.begin(), markerI.end(), 0);

    return runGibbs(analysis, std::move(markerI));
}

int Analysis::runGibbs(AnalysisGraph *analysis, const std::vector<unsigned int> &markers)
{
    const auto copy(markers);
    return runGibbs(analysis, std::move(copy));
}

bool Analysis::didAccumulate() const
{
#if defined(MPI_ENABLED)
    if (!m_opt->useHybridMpi)
        return false;

    int didAccumulate = false;
    const int localDidAccumulate = m_didAccumulate;
    MPI_Allreduce(&localDidAccumulate, &didAccumulate, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD);
    return didAccumulate;
#else
    return false;
#endif
}

void Analysis::accumulate(const KernelPtr &kernel, const ConstAsyncResultPtr &result)
{
    (void) kernel; // Unused
    (void) result; // Unused
    m_didAccumulate = true;
}

void Analysis::resetAccumulators()
{
    m_didAccumulate = false;
}
