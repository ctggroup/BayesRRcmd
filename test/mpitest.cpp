
#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <iostream>
#include <iterator>
#include <vector>

#include <thread>

int main(int argc, const char ** argv) {
#ifdef MPI_ENABLED

//    MPI_Init(nullptr, nullptr);

    int provided = 0;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);

    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << "MPI support enabled at compile time." << std::endl
                  << "World size: " << worldSize << std::endl;
    }

    std::cout << "Rank " << rank << " requested MPI thread support: "
              << MPI_THREAD_MULTIPLE << "; received: " << provided << std::endl;

    const long localRank = rank;
    std::vector<long> worldSizes(static_cast<size_t>(worldSize), 0);

    if (rank == 0) std::cout << "Beginning MPI_Barrier... " << std::endl;
    const auto barrierResult = MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::cout << " - complete! (" << barrierResult << ")" << std::endl;

    if (rank == 0) std::cout << "Beginning MPI_Allgather... " << std::endl;
    const auto allGatherResult = MPI_Allgather(&localRank, 1, MPI_LONG, worldSizes.data(), 1, MPI_LONG, MPI_COMM_WORLD);
    if (rank == 0) std::cout << " - complete! (" << allGatherResult << ")" << std::endl;

    std::cout << rank << ": ";
    std::copy(worldSizes.begin(), std::prev(worldSizes.end()), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << worldSizes.back() << std::endl;

    MPI_Finalize();
    return 0;
#else
    std::cout << "MPI support was not available at compile time." << std::endl;
    return 1;
#endif
}
