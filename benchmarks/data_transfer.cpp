#include <chrono>
#include <Eigen/Eigen>
#include <iostream>
#include <iterator>
#include <mpi.h>

using namespace Eigen;

struct BenchOptions
{
    size_t iterations = 100;
    long size = 500000;

    bool parse(int argc, const char* argv[]) {

        for (int i = 1; i < argc; ++i) {
            if (!strcmp(argv[i], "--size")) {
                size = std::stol(argv[++i]);
            } else if (!strcmp(argv[i], "--iterations")) {
                iterations = std::stoul(argv[++i]);
            }
        }

        return true;
    }
};

void benchmark(const BenchOptions& options, VectorXd& data) {
    VectorXd destination(options.size);
    MPI_Allreduce(data.data(), destination.data(), options.size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

int main(int argc, const char* argv[]) {
    BenchOptions options;
    options.parse(argc, argv);

    MPI_Init(nullptr, nullptr);

    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        std::cout << "Running benchmark with:" << std::endl
                  << "iterations: " << options.iterations << std::endl
                  << "size: " << options.size << " doubles" << std::endl
                  << "hosts: " << worldSize << std::endl;
    }

    std::chrono::milliseconds duration(0);
    VectorXd data = VectorXd::Ones(options.size);\

    for (size_t i = 0; i < options.iterations; ++i) {
        const auto start = std::chrono::steady_clock::now();

        benchmark(options, data);

        const auto end = std::chrono::steady_clock::now();
        duration += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    }

    std::chrono::milliseconds meanDuration = duration / options.iterations;
    const auto count = meanDuration.count();

    std::vector<long> results(worldSize);

    MPI_Gather(&count, 1, MPI_LONG, results.data(), 1, MPI_LONG, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Mean benchmark duration (ms): ";
        std::copy(results.begin(), std::prev(results.end()), std::ostream_iterator<long>(std::cout, ", "));
        std::cout << results.back() << std::endl;
    }

    MPI_Finalize();
    return 0;
}
