#ifndef HYBRIDMPISYNCMANAGER_H
#define HYBRIDMPISYNCMANAGER_H

#include "common.h"

#include <shared_mutex>

using SynchronisedMarkers = std::vector<unsigned int>;

namespace HybridMpi {

    using IndexResultPair = std::pair<unsigned int, AsyncResultPtr>;
    using ConstIndexResultPair = std::pair<unsigned int, ConstAsyncResultPtr>;

    class SyncManager
    {
    public:
        explicit SyncManager(bool useMpi);

        void accumulate(const KernelPtr& kernel, const ConstAsyncResultPtr& result);
        SynchronisedMarkers synchroniseMarkers();
        AsyncResultPtr resultForMarker(unsigned int marker) const;
        void reset();

    private:
        bool m_useMpi = true;
        int m_rank = 0;
        int m_worldSize = 1;

        mutable std::shared_mutex m_accumulatorMutex;
        std::vector<ConstIndexResultPair> m_localResults;

        std::vector<IndexResultPair> m_remoteResults;

#if defined(MPI_ENABLED) && defined(MPI_TIMING_ENABLED)
        double m_accumulateTime;
#endif
    };

    std::streamsize sendSize(const std::vector<ConstIndexResultPair>& pairs);
    std::vector<IndexResultPair> readResultPairs(std::istream *inStream);
    void writeResultPairs(const std::vector<ConstIndexResultPair>& pairs, std::ostream *outStream);

    std::streamsize size(const ConstAsyncResultPtr& result);
    void read(AsyncResultPtr& result, std::istream *inStream);
    void write(const ConstAsyncResultPtr& result, std::ostream *outStream);

    bool operator==(const AsyncResult &lhs, const AsyncResult &rhs);

}

#endif // HYBRIDMPISYNCMANAGER_H
