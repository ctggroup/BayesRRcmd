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

void processDenseData(Options opt) {
    Data data;

    const auto ppFile = ppFileForType(opt.dataType, opt.dataFile);
    const auto ppIndexFile = ppIndexFileForType(opt.dataType, opt.dataFile);
    // Analysis type
    if (opt.analysisType == "RAMBayes" && ( opt.bayesType == "bayes" || opt.bayesType == "bayesMmap" || opt.bayesType == "horseshoe"))
    {
        //TODO reimplement RAM and other solutions in a cleaner manner
        throw( "Error " + opt.analysisType +  " not implemented yet");
    }
    else if(opt.analysisType == "Preprocess") //preprocess and compress data
    {
        cout << "Start preprocessing " << opt.dataFile << endl;
        clock_t start_bed = clock();

        if(opt.inputType == InputType::BED)
        {
            data.readFamFile(fileWithSuffix(opt.dataFile, ".fam"));
            data.readBimFile(fileWithSuffix(opt.dataFile, ".bim"));
            std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
            if (opt.numThreadSpawned > 0)
                taskScheduler = std::make_unique<tbb::task_scheduler_init>(opt.numThreadSpawned);

            std::cout << "Preprocessing with " << opt.numThread << " threads ("
                      << (opt.numThreadSpawned > 0 ? std::to_string(opt.numThreadSpawned) : "auto") << " spawned) and "
                      << opt.preprocessChunks << " columns per thread."
                      << endl;

            PreprocessGraph graph(opt.numThread);
            graph.preprocessBedFile(opt.dataFile,
                                    opt.dataType,
                                    opt.compress,
                                    &data,
                                    opt.preprocessChunks);
        }
        else
        {
            data.readCSVFile(opt.dataFile);
            data.preprocessCSVFile(opt.dataFile, ppFile, ppIndexFile, opt.compress);
        } //end if bedfile

        clock_t end = clock();
        printf("Finished preprocessing in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
        cout << endl;

    }//end if preprocess
    else if(opt.analysisType == "PPBayes" || opt.analysisType == "PPAsyncBayes")
    {
        if(opt.inputType == InputType::BED)
        {
            data.readFamFile(fileWithSuffix(opt.dataFile, ".fam"));
            data.readBimFile(fileWithSuffix(opt.dataFile, ".bim"));
            data.readPhenotypeFile(opt.phenotypeFile);
        }
        else
        {
            data.readCSVFile(opt.dataFile);
            data.readCSVPhenFile(opt.phenotypeFile);
        }

        // Run analysis using mapped data files
        cout << "Start reading preprocessed bed file: " << ppFile << endl;
        clock_t start_bed = clock();
        data.mapCompressedPreprocessBedFile(ppFile, ppIndexFile);
        clock_t end = clock();
        printf("Finished reading preprocessed bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
        cout << endl;

        std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
        if (opt.numThreadSpawned > 0)
            taskScheduler = std::make_unique<tbb::task_scheduler_init>(opt.numThreadSpawned);

        std::unique_ptr<AnalysisGraph> graph {nullptr};
        if (opt.analysisType == "PPAsyncBayes")
        {
            graph = std::make_unique<ParallelGraph>(opt.decompressionTokens, opt.analysisTokens);
            auto *parallelGraph = dynamic_cast<ParallelGraph*>(graph.get());
            parallelGraph->setDecompressionNodeConcurrency(opt.decompressionNodeConcurrency);
            parallelGraph->setAnalysisNodeConcurrency(opt.analysisNodeConcurrency);
        }
        else
        {
            graph = std::make_unique<LimitSequenceGraph>(opt.numThread);
        }
        DenseBayesRRmz analysis(&data, opt);
        analysis.runGibbs(graph.get());
        data.unmapCompressedPreprocessedBedFile();
    }//end if ppbayes
    else
    {
        throw(" Error: Wrong analysis type: " + opt.analysisType);
    }//end if analysis type
}

void processSparseData(Options options) {
    if (options.analysisType != "PPBayes" &&
            options.analysisType != "PPAsyncBayes" &&
            options.analysisType != "Preprocess") {
        std::cout << "Error: Wrong analysis type: " << options.analysisType << std::endl;
        return;
    }

    Data data;

    // Read in the data for every possible option

    data.readFamFile(fileWithSuffix(options.dataFile, ".fam"));
    data.readBimFile(fileWithSuffix(options.dataFile, ".bim"));
    data.readPhenotypeFile(options.phenotypeFile);

    const auto ppFile = ppFileForType(options.dataType, options.dataFile);
    const auto ppIndexFile = ppIndexFileForType(options.dataType, options.dataFile);

    if (options.analysisType == "Preprocess") {
        cout << "Start preprocessing " << options.dataFile << endl;

        clock_t start_bed = clock();

        std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
        if (options.numThreadSpawned > 0)
            taskScheduler = std::make_unique<tbb::task_scheduler_init>(options.numThreadSpawned);

        std::cout << "Preprocessing with " << options.numThread << " threads ("
                  << (options.numThreadSpawned > 0 ? std::to_string(options.numThreadSpawned) : "auto") << " spawned) and "
                  << options.preprocessChunks << " columns per thread."
                  << endl;


        PreprocessGraph graph(options.numThread);
        graph.preprocessBedFile(options.dataFile,
                                options.dataType,
                                options.compress,
                                &data,
                                options.preprocessChunks);

        clock_t end = clock();
        printf("Finished preprocessing the bed file in %.3f sec.\n",
               double(end - start_bed) / double(CLOCKS_PER_SEC));
        cout << endl;
        return;
    }

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
    if (options.analysisType == "PPAsyncBayes") {
        graph = std::make_unique<ParallelGraph>(options.decompressionTokens, options.analysisTokens);
        auto *parallelGraph = dynamic_cast<ParallelGraph*>(graph.get());
        parallelGraph->setDecompressionNodeConcurrency(options.decompressionNodeConcurrency);
        parallelGraph->setAnalysisNodeConcurrency(options.analysisNodeConcurrency);
    } else {
        graph = std::make_unique<LimitSequenceGraph>(options.numThread);
    }

    SparseBayesRRG analysis(&data, options);
    analysis.runGibbs(graph.get());

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
        Options opt;
        opt.inputOptions(argc, argv);

        if (opt.inputType == InputType::Unknown) {
            cerr << "Unknown --datafile type: '" << opt.dataFile << "'\n";
            exit(1);
        }

        switch (opt.dataType) {
        case DataType::Dense:
            processDenseData(opt);
            break;

        case DataType::SparseEigen:
            // Fall through
        case DataType::SparseRagged:
            processSparseData(opt);
            break;

        default:
            cerr << "Unsupported DataType: " << opt.dataType << endl;
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
