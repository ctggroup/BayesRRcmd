#include <iostream>
#include <string>
#include "BayesRRm.h"
#include "DenseBayesRRmz.hpp"
#include "data.hpp"
#include "options.hpp"
#include "SparseBayesRRG.hpp"
#include "tbb/task_scheduler_init.h"
#include "preprocessgraph.h"
#include "densemarker.h"
#include "raggedsparsemarker.h"
#include "common.h"
#include "limitsequencegraph.hpp"
#include "parallelgraph.h"

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
};

void preprocessBed(const Options &options) {
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
}

void preprocessCsv(const Options &options) {
    assert(options.analysisType == AnalysisType::Preprocess);
    assert(options.inputType == InputType::CSV);

    if (options.preprocessDataType != PreprocessDataType::Dense) {
        cerr << "CSV preprocessing only supports the Dense data type." << endl;
        return;
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
}

void preprocess(const Options &options) {
    assert(options.analysisType == AnalysisType::Preprocess);

    switch (options.inputType) {
    case InputType::BED:
        preprocessBed(options);
        break;

    case InputType::CSV:
        preprocessCsv(options);
        break;

    default:
        cout << "Cannot preprocess for input type: " << options.inputType << endl;
        return;
    }
}

void runPpBayesAnalysis(const Options &options) {
    assert(options.analysisType == AnalysisType::PpBayes ||
           options.analysisType == AnalysisType::AsyncPpBayes);

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

    std::unique_ptr<AnalysisGraph> graph {nullptr};
    if (options.analysisType == AnalysisType::AsyncPpBayes) {
        graph = std::make_unique<ParallelGraph>(options.decompressionTokens, options.analysisTokens);
        auto *parallelGraph = dynamic_cast<ParallelGraph*>(graph.get());
        parallelGraph->setDecompressionNodeConcurrency(options.decompressionNodeConcurrency);
        parallelGraph->setAnalysisNodeConcurrency(options.analysisNodeConcurrency);
    } else {
        graph = std::make_unique<LimitSequenceGraph>(options.numThread);
    }

    switch (options.preprocessDataType) {
    case PreprocessDataType::Dense:
    {
        DenseBayesRRmz analysis(&data, options);
        analysis.runGibbs(graph.get());
        break;
    }

    case PreprocessDataType::SparseEigen:
        // Fall through
    case PreprocessDataType::SparseRagged:
    {
        SparseBayesRRG analysis(&data, options);
        analysis.runGibbs(graph.get());
        break;
    }

    default:
        cerr << "Unsupported DataType: " << options.preprocessDataType << endl;
        break;
    }

    data.unmapCompressedPreprocessedBedFile();
}

int main(int argc, const char * argv[]) {


    cout << "***********************************************\n";
    cout << "* BayesRRcmd                                  *\n";
    cout << "* Complex Trait Genetics group UNIL           *\n";
    cout << "*                                             *\n";
    cout << "* MIT License                                 *\n";
    cout << "***********************************************\n";

    Gadget::Timer timer;
    timer.setTime();
    cout << "\nAnalysis started: " << timer.getDate();

    if (argc < 2){
        cerr << " \nDid you forget to give the input parameters?\n" << endl;
        exit(1);
    }
    try {
        Options options;
        options.inputOptions(argc, argv);

        if (options.inputType == InputType::Unknown) {
            cerr << "Unknown --datafile type: '" << options.dataFile << "'\n";
            exit(1);
        }

        switch (options.analysisType) {
        case AnalysisType::Preprocess:
            preprocess(options);
            break;

        case AnalysisType::PpBayes:
            // Fall through
        case AnalysisType::AsyncPpBayes:
            runPpBayesAnalysis(options);
            break;

        default:
            cerr << "Unknown --analyis-type" << endl;
            exit(1);
        }
    }
    catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    }
    catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }

    timer.getTime();

    cout << "\nAnalysis finished: " << timer.getDate();
    cout << "Computational time: "  << timer.format(timer.getElapse()) << endl;


    return 0;
}
