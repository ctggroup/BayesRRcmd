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

using namespace std;

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
    readMetaData(data, options);

    PreprocessGraph graph(options.numThread);
    graph.preprocessBedFile(options.dataFile,
                            options.preprocessDataType,
                            options.compress,
                            &data,
                            options.preprocessChunks);

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
    readMetaData(data, options);

    const auto ppFile = ppFileForType(options.preprocessDataType, options.dataFile);
    const auto ppIndexFile = ppIndexFileForType(options.preprocessDataType, options.dataFile);
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
switch (options->preprocessDataType) {
    case PreprocessDataType::Dense:
    {
        DenseBayesRRmz analysis(data, options);
        analysis.runGibbs(graph);
        break;
    }

    case PreprocessDataType::SparseEigen:
        // Fall through
    case PreprocessDataType::SparseRagged:
    {
        SparseBayesRRG analysis(data, options);
        analysis.runGibbs(graph);
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

    switch (options->preprocessDataType) {
    case PreprocessDataType::Dense:
    {
        DenseBayesW analysis(data, options, sysconf(_SC_PAGE_SIZE));
        analysis.runGibbs(graph);
        break;
    }

    case PreprocessDataType::SparseRagged:
    {
        SparseBayesW analysis(data, options, sysconf(_SC_PAGE_SIZE));
        analysis.runGibbs(graph);
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
    readMetaData(data, options);

    const auto ppFile = ppFileForType(options.preprocessDataType, options.dataFile);
    const auto ppIndexFile = ppIndexFileForType(options.preprocessDataType, options.dataFile);

    cout << "Start reading preprocessed bed file: " << ppFile << endl;
    clock_t start_bed = clock();
    data.mapCompressedPreprocessBedFile(ppFile, ppIndexFile);
    clock_t end = clock();
    printf("Finished reading preprocessed bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
    cout << endl;

    std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
    if (options.numThreadSpawned > 0)
        taskScheduler = std::make_unique<tbb::task_scheduler_init>(options.numThreadSpawned);

    if (options.useMarkerCache) {
        clock_t start_cache = clock();
        markerCache()->populate(&data, &options);
        clock_t end_cache = clock();
        printf("Populated cache in %.3f sec.\n", double(end_cache - start_cache) / double(CLOCKS_PER_SEC));
        cout << endl;
    }

    auto graph = AnalysisRunner::makeAnalysisGraph(options);

    auto cleanup = [&data]() {
        data.unmapCompressedPreprocessedBedFile();
    };

    bool result = false;
    switch (options.analysisType) {
    case AnalysisType::PpBayes:
        // Fall through
    case AnalysisType::AsyncPpBayes:
        result = runBayesRAnalysis(&options, &data, graph.get());
        break;

    case AnalysisType::Gauss:
        // Fall through
    case AnalysisType::AsyncGauss:
        result = runBayesWAnalysis(&options, &data, graph.get());
        break;

    default:
        cerr << "Unsupported AnalysisType: runBayesAnalysis does not support"
             << options.analysisType << endl;
        cleanup();
        return false;
    }

    cleanup();
    return result;
}

namespace AnalysisRunner {

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
                                                             options.useMarkerCache);
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

    default:
        cerr << "Unknown --analyis-type" << endl;
        return false;
    }
}

}
