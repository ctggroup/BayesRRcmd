#include "analysisrunner.h"

#include "tbb/task_scheduler_init.h"

#include "common.h"
#include "data.hpp"
#include "DenseBayesRRmz.hpp"
#include "densebayesw.h"
#include "limitsequencegraph.hpp"
#include "markercache.h"
#include "options.hpp"
#include "parallelgraph.h"
#include "preprocessgraph.h"
#include "SparseBayesRRG.hpp"
#include "sequential.h"
#include "sparsebayesw.h"
#include "preprocessedfilesplitter.h"

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <utility>

using namespace std;

namespace {

class MpiContext {
public:
    MpiContext(bool useMpi)
        : m_useMpi(useMpi)
    {
#ifdef MPI_ENABLED
        if (m_useMpi) {
            MPI_Init(nullptr, nullptr);
            MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI::ERRORS_THROW_EXCEPTIONS);
        }
#endif
    }

    ~MpiContext() {
#ifdef MPI_ENABLED
        if (m_useMpi)
            MPI_Finalize();
#endif
    }

private:
    bool m_useMpi = false;
};

MarkerSubset getMarkerSubset(const Options *options, const Data *data) {
#ifdef MPI_ENABLED
    if (options->useHybridMpi) {
        int worldSize;
        MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

        const auto subsets = generateEqualSubsets(worldSize, data->numSnps);

        if (!isValid(subsets, data->numSnps)) {
            cerr << "Invalid marker distribution" << endl;
            for (const auto &s : subsets)
                cerr << s.first() << ", " << s.last() << ": " << s.size() << endl;
            MPI_Abort(MPI_COMM_WORLD, -2);
        }

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const auto subset = subsets.at(rank);
        if (subset.isValid(data->numSnps)) {
            printf("Rank %d working markers %d to %lo\n", rank, subset.first(), subset.last());
        } else {
            cerr << "Rank " << rank << " has invalid marker subset: " << subset.first() << " to " << subset.last()
                 << " for " << data->numSnps << " markers" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        return subset;
    } else
#endif
    {
        return options->preprocessSubset;
    }
}

bool split(const Options& options) {
    assert(options.analysisType == AnalysisType::Split);

    Data data;
    AnalysisRunner::readMetaData(data, options);
    data.mapPreprocessBedFile(options);

    PreprocessedFileSplitter splitter;
    return splitter.split(options, &data, getMarkerSubset(&options, &data));
}

bool preprocessBed(const Options &options) {
    assert(options.analysisType == AnalysisType::Preprocess);
    assert(options.inputType == InputType::BED);

    cout << "Start preprocessing " << options.dataFile << endl
         << "Preprocessing with " << options.numThread << " threads ("
         << (options.numThreadSpawned > 0 ? std::to_string(options.numThreadSpawned) : "auto") << " spawned) and "
         << options.preprocessChunks << " columns per thread."
         << endl;

    clock_t start_bed = clock();

    std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
    if (options.numThreadSpawned > 0)
        taskScheduler = std::make_unique<tbb::task_scheduler_init>(options.numThreadSpawned);

    Data data;
    AnalysisRunner::readMetaData(data, options);

    data.setMarkerSubset(getMarkerSubset(&options, &data));

    PreprocessGraph graph(options.numThread);
    graph.preprocessBedFile(options, &data);

    clock_t end = clock();
    printf("Finished preprocessing the bed file in %.3f sec.\n\n",
           double(end - start_bed) / double(CLOCKS_PER_SEC));

    return true;
}

bool preprocessCsv(const Options &options) {
    assert(options.analysisType == AnalysisType::Preprocess);
    assert(options.inputType == InputType::CSV);

    if (options.preprocessDataType != PreprocessDataType::Dense) {
        cerr << "CSV preprocessing only supports the Dense data type." << endl;
        return false;
    }

    cout << "Start preprocessing " << options.dataFile << endl;

    clock_t start_bed = clock();

    Data data;
    AnalysisRunner::readMetaData(data, options);

    const auto ppFile = ppFileForType(options);
    const auto ppIndexFile = ppIndexFileForType(options);
    data.preprocessCSVFile(options.dataFile, ppFile, ppIndexFile, options.compress);

    clock_t end = clock();
    printf("Finished preprocessing the bed file in %.3f sec.\n\n",
           double(end - start_bed) / double(CLOCKS_PER_SEC));

    return true;
}

bool preprocess(const Options &options) {
    assert(options.analysisType == AnalysisType::Preprocess);

    switch (options.inputType) {
    case InputType::BED:
        return preprocessBed(options);

    case InputType::CSV:
        return preprocessCsv(options);

    default:
        cout << "Cannot preprocess for input type: " << options.inputType << endl;
        return false;
    }
}

bool runBayesRAnalysis(const Options *options, Data *data, AnalysisGraph *graph) {
    // If there is a file for fixed effects (model matrix), then read the data
    if(!options->fixedFile.empty()) {
        data->readCSV(options->fixedFile, options->fixedEffectNumber);
    }

    auto markers = data->getMarkerIndexList();

    switch (options->preprocessDataType) {
    case PreprocessDataType::Dense:
    {
        DenseBayesRRmz analysis(data, options);
        analysis.runGibbs(graph, std::move(markers));
        break;
    }

    case PreprocessDataType::SparseEigen:
        // Fall through
    case PreprocessDataType::SparseRagged:
    {
        SparseBayesRRG analysis(data, options);
        analysis.runGibbs(graph, std::move(markers));
        break;
    }

    default:
        cerr << "Unsupported DataType: BayesR does not support"
             << options->preprocessDataType << endl;
        return false;
    }

    return true;
}

bool runBayesWAnalysis(const Options *options, Data *data, AnalysisGraph *graph) {
    // If there is a file for fixed effects (model matrix), then read the data
    if(!options->fixedFile.empty()) {
        data->readCSV(options->fixedFile, options->fixedEffectNumber);
    }

    // Read the failure indicator vector
    data->readFailureFile(options->failureFile);

    auto markers = data->getMarkerIndexList();

    switch (options->preprocessDataType) {
    case PreprocessDataType::Dense:
    {
        DenseBayesW analysis(data, options, sysconf(_SC_PAGE_SIZE));
        analysis.runGibbs(graph, std::move(markers));
        break;
    }

    case PreprocessDataType::SparseRagged:
    {
        SparseBayesW analysis(data, options, sysconf(_SC_PAGE_SIZE));
        analysis.runGibbs(graph, std::move(markers));
        break;
    }

    default:
        cerr << "Unsupported DataType: BayesW does not support"
             << options->preprocessDataType << endl;
        return false;
    }

    return true;
}


bool runBayesAnalysis(const Options &options) {
    Data data;
    AnalysisRunner::readMetaData(data, options);

    cout << "Start reading preprocessed bed file: " << ppFileForType(options) << endl;
    clock_t start_bed = clock();
    data.mapPreprocessBedFile(options);
    clock_t end = clock();
    printf("Finished reading preprocessed bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
    cout << endl;

    if (!data.validMarkerSubset()) {
        const auto first = data.markerSubset().first();
        const auto last = data.markerSubset().last();
        cerr << "Marker range " << first  << " to " << last << "is not valid!" << endl
             << "Expected range is 0 to " << data.numSnps - 1 << endl;
        return false;
    }

    std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
    if (options.numThreadSpawned > 0)
        taskScheduler = std::make_unique<tbb::task_scheduler_init>(options.numThreadSpawned);

    if (options.useMarkerCache) {
        markerCache()->populate(&data, &options);
    }

    auto graph = AnalysisRunner::makeAnalysisGraph(options);

    switch (options.analysisType) {
    case AnalysisType::PpBayes:
        // Fall through
    case AnalysisType::AsyncPpBayes:
        return runBayesRAnalysis(&options, &data, graph.get());
        break;

    case AnalysisType::Gauss:
        // Fall through
    case AnalysisType::AsyncGauss:
        return runBayesWAnalysis(&options, &data, graph.get());
        break;

    default:
        cerr << "Unsupported AnalysisType: runBayesAnalysis does not support"
             << options.analysisType << endl;
        return false;
    }
}

}

namespace AnalysisRunner {

void readMetaData(Data &data, const Options &options) {
    switch (options.inputType) {
    case InputType::BED:
        data.readFamFile(fileWithSuffix(options.dataFile, ".fam"));
        data.readBimFile(fileWithSuffix(options.dataFile, ".bim"));
        data.readPhenotypeFile(options.phenotypeFile);
        break;

    case InputType::CSV:
        data.readCSVFile(options.dataFile);
        data.readCSVPhenFile(options.phenotypeFile);
        break;

    default:
        cout << "Cannot read meta data for input type: " << options.inputType << endl;
        return;
    }

    if (!options.groupFile.empty()) {
        data.readGroupFile(options.groupFile);
        if (options.S.rows() != data.numGroups)
            cerr << "Number of groups " << data.numGroups
                 << " does not match the number of variance sets: " << options.S.rows()
                 << endl;
    }
};

std::unique_ptr<AnalysisGraph> makeAnalysisGraph(const Options &options)
{
    switch (options.analysisType) {
    case AnalysisType::PpBayes:
        // Fall through
    case AnalysisType::Gauss:
    {
        if (options.useMarkerCache)
            return std::make_unique<::Sequential>(); // Differentiate from Eigen::Sequential
        else
            return std::make_unique<LimitSequenceGraph>(options.numThread);
    }

    case AnalysisType::AsyncPpBayes:
        // Fall through
    case AnalysisType::AsyncGauss:
    {
        auto parallelGraph = std::make_unique<ParallelGraph>(options.decompressionTokens,
                                                             options.analysisTokens,
                                                             options.useMarkerCache,
                                                             options.useHybridMpi);
        parallelGraph->setDecompressionNodeConcurrency(options.decompressionNodeConcurrency);
        parallelGraph->setAnalysisNodeConcurrency(options.analysisNodeConcurrency);
        return std::move(parallelGraph);
    }

    default:
        std::cerr << "makeAnalysisGraph - unsupported AnalysisType: " << options.analysisType
                  << std::endl;
        return {};
    }
}

bool run(const Options &options)
{
    if (options.inputType == InputType::Unknown) {
        cerr << "Unknown --datafile type: '" << options.dataFile << "'" << endl;
        return false;
    }

    if (!options.validWorkingDirectory()) {
        return false;
    }

#if !defined(MPI_ENABLED)
    if (options.useHybridMpi) {
        cerr << "--hybrid-mpi flag has been used but MPI support was not compiled in!" << endl;
        return false;
    }
#endif

    MpiContext mpiContext(options.useHybridMpi);

    bool result = false;
    try {
        switch (options.analysisType) {
        case AnalysisType::Preprocess:
            return preprocess(options);

        case AnalysisType::PpBayes:
            // Fall through
        case AnalysisType::AsyncPpBayes:
            // Fall through
        case AnalysisType::Gauss:
            // Fall through
        case AnalysisType::AsyncGauss:
            if (options.S.size() == 0) {
                cerr << "Variance components `--S` are missing or could not be parsed!" << endl;
                return false;
            }
            return runBayesAnalysis(options);

        case AnalysisType::Split:
            return split(options);

        default:
            cerr << "Unknown --analyis-type" << endl;
            break;
        }
    }
    catch (const string &err_msg) {
        cerr << endl << err_msg << endl;
#ifdef MPI_ENABLED
        if (options.useHybridMpi)
            MPI::COMM_WORLD.Abort(-1);
#endif
    }
    catch (const char *err_msg) {
        cerr << endl << err_msg << endl;
#ifdef MPI_ENABLED
        if (options.useHybridMpi)
            MPI::COMM_WORLD.Abort(-1);
#endif
    }
    catch (const std::exception &e) {
        cerr << endl << e.what() << endl;
#ifdef MPI_ENABLED
        if (options.useHybridMpi)
            MPI::COMM_WORLD.Abort(-1);
#endif
    }
#ifdef MPI_ENABLED
    catch (MPI::Exception e) {
        // Not all MPI implementations mark these functions as const
        cerr << endl << e.Get_error_string() << endl;
        if (options.useHybridMpi)
            MPI::COMM_WORLD.Abort(e.Get_error_code());
    }
#endif

    return false;
}

}
