#include "hybridmpisyncmanager.h"

#include "asyncresult.h"
#include "kernel.h"

#include <numeric>

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

namespace HybridMpi {

namespace  {

struct IndexResultCompareAdapter {
    IndexResultCompareAdapter(const ConstIndexResultPair & p) : i(p.first) {}
    IndexResultCompareAdapter(const IndexResultPair& p) : i(p.first) {}

    unsigned int i = 0;
};

bool operator<(const IndexResultCompareAdapter& lhs, const IndexResultCompareAdapter& rhs) {
    return lhs.i < rhs.i;
}

}

SyncManager::SyncManager(bool useMpi)
    : m_useMpi(useMpi)
{
#ifdef MPI_ENABLED
    if (m_useMpi) {
        MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &m_worldSize);
    }
#endif
}

void SyncManager::accumulate(const KernelPtr &kernel, const ConstAsyncResultPtr &result)
{
#ifdef MPI_TIMING_ENABLED
    const auto start = std::chrono::steady_clock::now();
#endif
    std::unique_lock lock(m_accumulatorMutex);
    m_localResults.emplace_back(kernel->marker->i, result);
#ifdef MPI_TIMING_ENABLED
    const auto end = std::chrono::steady_clock::now();
    m_accumulateTime += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000000.0;
#endif
}

SynchronisedMarkers SyncManager::synchroniseMarkers()
{
#if defined(MPI_ENABLED)
    if (!m_useMpi)
        return {};

    // Gather the sizes of data to transfer from each host
    const long localSize = sendSize(m_localResults);
    std::vector<long> worldSizes(m_worldSize, 0);

    MPI_Allgather(&localSize, 1, MPI_LONG, worldSizes.data(), 1, MPI_LONG, MPI_COMM_WORLD);

    const std::streamsize worldDataSize = std::accumulate(worldSizes.begin(), worldSizes.end(), 0);
    if (worldDataSize == 0)
        return {};

    // Sort our results in preparation for sending
    auto indexResultCompare = [](const IndexResultCompareAdapter& a, const IndexResultCompareAdapter& b) {
        return a < b;
    };
    std::sort(m_localResults.begin(), m_localResults.end(), indexResultCompare);

    // Create our send buffer
    std::unique_ptr<char[]> outBuffer;
    if (localSize > 0) {
        outBuffer.reset(new char[static_cast<unsigned long>(localSize)]);

        std::ostringstream stream;
        stream.rdbuf()->pubsetbuf(outBuffer.get(), localSize);

        writeResultPairs(m_localResults, &stream);
    }

    // Create our receive buffer, recvCount and displacement data
    std::unique_ptr<char[]> inBuffer;
    inBuffer.reset(new char[static_cast<unsigned long>(worldDataSize)]);

    std::vector<int> recvCount;
    std::transform(worldSizes.begin(), worldSizes.end(), std::back_inserter(recvCount), [](const long l) { return static_cast<int>(l); });

    std::vector<int> displacements = recvCount;
    std::rotate(displacements.rbegin(), displacements.rbegin() + 1, displacements.rend());
    displacements.at(0) = 0;

    // MPI_Allgather the results data between hosts
    MPI_Allgatherv(outBuffer.get(), localSize, MPI_CHAR, inBuffer.get(), recvCount.data(), displacements.data(), MPI_CHAR, MPI_COMM_WORLD);

    // Convert receive buffer to std::pair<unsigned int, AsyncResultPtr>
    std::vector<IndexResultPair> allResults;
    {
        std::istringstream stream;
        stream.rdbuf()->pubsetbuf(inBuffer.get(), worldDataSize);
        allResults = readResultPairs(&stream);
    }

    if (allResults.empty())
        return {};

    // Separate local results from all results
    std::sort(allResults.begin(), allResults.end(), indexResultCompare);
    m_remoteResults.clear();
    m_remoteResults.reserve(allResults.size());
    std::set_difference(allResults.begin(), allResults.end(),
                        m_localResults.begin(), m_localResults.end(),
                        std::inserter(m_remoteResults, m_remoteResults.begin()),
                        indexResultCompare);

    // Return the indexes of the markers that need processing
    SynchronisedMarkers markers;
    markers.reserve(m_remoteResults.size());
    std::transform(m_remoteResults.begin(), m_remoteResults.end(), std::back_inserter(markers), [](const IndexResultPair& p) {
        return p.first;
    });

    return markers;
#else
    return {};
#endif
}

void SyncManager::reset()
{
    m_localResults.clear();
}

std::streamsize sendSize(const std::vector<ConstIndexResultPair> &pairs)
{
    if (pairs.empty())
        return 0;

    const auto resultSize = size(pairs.front().second);
    return static_cast<std::streamsize>((sizeof (unsigned int) + resultSize) * pairs.size());
}

std::vector<IndexResultPair> readResultPairs(std::istream *inStream)
{
    std::vector<IndexResultPair> results;
    while (true) {
        unsigned int index = 0;
        inStream->read(reinterpret_cast<char *>(&index), sizeof (unsigned int));
        if (!inStream->good())
            break;

        AsyncResultPtr result;
        HybridMpi::read(result, inStream);
        if (!inStream->good())
            break;

        results.emplace_back(index, result);
    }

    if (!inStream->eof()) {
        std::cerr << "Failed to read IndexResultPairs for HybridMpi::SyncManager!" << std::endl;
        return {};
    }

    return results;
}

void writeResultPairs(const std::vector<ConstIndexResultPair> &pairs, std::ostream *outStream)
{
    std::for_each(pairs.begin(), pairs.end(), [&outStream](const ConstIndexResultPair& pair) {
        outStream->write(reinterpret_cast<const char *>(&pair.first), sizeof (unsigned int));
        HybridMpi::write(pair.second, outStream);
    });
}

std::streamsize size(const ConstAsyncResultPtr &result)
{
    return static_cast<std::streamsize>(sizeof (result->betaOld)
                                        + sizeof (result->beta)
                                        + sizeof (result->v->rows())
                                        + sizeof (result->v->cols())
                                        + (result->v->size() * sizeof (double)) // v data
                                        + sizeof (result->component));
}

void read(AsyncResultPtr &result, std::istream *inStream)
{
    if (inStream->fail()) {
        std::cerr << "Error: unable to read AsyncResult for HybridMpi::SyncManager!" << std::endl;
        return;
    }

    result.reset(new AsyncResult);

    auto readDouble = [&inStream](double *d) {
        inStream->read(reinterpret_cast<char *>(d), sizeof(double));
    };

    auto readIndex = [&inStream](Eigen::Index *i) {
        inStream->read(reinterpret_cast<char *>(i), sizeof(Eigen::Index));
    };

    readDouble(&result->betaOld);
    readDouble(&result->beta);

    Eigen::Index rows = 0;
    readIndex(&rows);

    Eigen::Index cols = 0;
    readIndex(&cols);

    result->v = std::make_unique<MatrixXd>(rows, cols);
    inStream->read(reinterpret_cast<char *>(result->v->data()),
                   static_cast<std::streamsize>(result->v->size() * sizeof (double)));

    readDouble(&result->component);
}

void write(const ConstAsyncResultPtr &result, std::ostream *outStream)
{
    if (outStream->fail()) {
        std::cerr << "Error: unable to write AsyncResult for HybridMpi::SyncManager!" << std::endl;
        return;
    }

    auto writeDouble = [&outStream](const double d) {
        outStream->write(reinterpret_cast<const char *>(&d), sizeof(double));
    };

    auto writeIndex = [&outStream](const Eigen::Index i) {
        outStream->write(reinterpret_cast<const char *>(&i), sizeof(Eigen::Index));
    };

    writeDouble(result->betaOld);
    writeDouble(result->beta);
    writeIndex(result->v->rows());
    writeIndex(result->v->cols());
    outStream->write(reinterpret_cast<const char *>(result->v->data()),
                     static_cast<std::streamsize>(result->v->size() * sizeof (double)));
    writeDouble(result->component);
}

auto asTuple(const AsyncResult &r) {
    return std::tie(r.betaOld, r.beta, *(r.v), r.component);
};

bool operator==(const AsyncResult &lhs, const AsyncResult &rhs)
{
    return asTuple(lhs) == asTuple(rhs);
}

}
