#include <iostream>
#include <string>
#include "BayesRRm.h"
#include "BayesW.hpp"
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

    // Read in the data for every possible option
    data.readFamFile(opt.bedFile + ".fam");
    data.readBimFile(opt.bedFile + ".bim");

    const auto bedFile = opt.bedFile + ".bed";
    const auto ppFile = ppFileForType(opt.dataType, opt.bedFile);
    const auto ppIndexFile = ppIndexFileForType(opt.dataType, opt.bedFile);

    // RAM solution (analysisType = RAMBayes)

    if (opt.analysisType == "RAMBayes" && ( opt.bayesType == "bayes" || opt.bayesType == "bayesMmap" || opt.bayesType == "horseshoe" || opt.bayesType == "bayesW")) {

        clock_t start = clock();

        // Read phenotype file and bed file for the option specified
        data.readPhenotypeFile(opt.phenotypeFile);
        data.readBedFile_noMPI(opt.bedFile+".bed");
        data.readBedFile_noMPI_unstandardised(opt.bedFile+".bed"); // This part to read the non-standardised data

        // If there is a file for fixed effects (model matrix), then read the data
        if(opt.fixedFile != ""){
        	data.readCSV(opt.fixedFile,opt.fixedEffectNumber);
        }

        // Option bayesType="bayesMmap" is going to be deprecated
        if (opt.bayesType == "bayesMmap" || opt.bayesType == "bayes"){
            BayesRRm analysis(data, opt, sysconf(_SC_PAGE_SIZE));
            analysis.runGibbs();
        } else if (opt.bayesType == "horseshoe") {
            //TODO Finish horseshoe
        } else if (opt.bayesType == "bayesW") {
            BayesW analysis(data, opt, sysconf(_SC_PAGE_SIZE));
            if(opt.bayesW_version == "gauss"){
            	analysis.runGibbs_Gauss();
            }else{
            	cout << "Choose either bayesW_version = gauss" << endl;
            }

        } else if (opt.bayesType == "bayesG") {
            //TODO add Bayes groups
        }

        clock_t end   = clock();
        printf("OVERALL read+compute time = %.3f sec.\n", (float)(end - start) / CLOCKS_PER_SEC);
    }

    // Pre-processing the data (centering and scaling)
    else if (opt.analysisType == "Preprocess") {
        cout << "Start preprocessing " << opt.bedFile + ".bed" << endl;

        clock_t start_bed = clock();
        if (opt.numThread > 1) {
            std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
            if (opt.numThreadSpawned > 0)
                taskScheduler = std::make_unique<tbb::task_scheduler_init>(opt.numThreadSpawned);

            std::cout << "Preprocessing with " << opt.numThread << " threads ("
                      << (opt.numThreadSpawned > 0 ? std::to_string(opt.numThreadSpawned) : "auto") << " spawned) and "
                      << opt.preprocessChunks << " columns per thread."
                      << endl;

            PreprocessGraph graph(opt.numThread);
            graph.preprocessBedFile(opt.bedFile,
                                    opt.dataType,
                                    opt.compress,
                                    &data,
                                    opt.preprocessChunks);
        } else {
            data.preprocessBedFile(bedFile, ppFile, ppIndexFile, opt.compress);
        }

        clock_t end = clock();
        printf("Finished preprocessing the bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
        cout << endl;
    }else if (opt.analysisType == "PPBayes" || opt.analysisType == "PPAsyncBayes") {
        clock_t start = clock();
        data.readPhenotypeFile(opt.phenotypeFile);
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
        if (opt.analysisType == "PPAsyncBayes") {
            graph = std::make_unique<ParallelGraph>(opt.numThread);
        } else {
            graph = std::make_unique<LimitSequenceGraph>(opt.numThread);
        }

        if(opt.bayesType == "bayesW"){
        	cout << "bayesW pp solution not completed" << endl;
        	exit(55);
        	//TODO
        }else{
            DenseBayesRRmz analysis(&data, opt);
            analysis.runGibbs(graph.get());
        }


        data.unmapCompressedPreprocessedBedFile();
    }else {
        throw(" Error: Wrong analysis type: " + opt.analysisType);
    }
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
    data.readFamFile(options.bedFile + ".fam");
    data.readBimFile(options.bedFile + ".bim");
    data.readPhenotypeFile(options.phenotypeFile);

    const auto bedFile = options.bedFile + ".bed";
    const auto ppFile = ppFileForType(options.dataType, options.bedFile);
    const auto ppIndexFile = ppIndexFileForType(options.dataType, options.bedFile);

    if (options.analysisType == "Preprocess") {
        cout << "Start preprocessing " << bedFile << endl;

        clock_t start_bed = clock();

        std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
        if (options.numThreadSpawned > 0)
            taskScheduler = std::make_unique<tbb::task_scheduler_init>(options.numThreadSpawned);

        std::cout << "Preprocessing with " << options.numThread << " threads ("
                  << (options.numThreadSpawned > 0 ? std::to_string(options.numThreadSpawned) : "auto") << " spawned) and "
                  << options.preprocessChunks << " columns per thread."
                  << endl;


        PreprocessGraph graph(options.numThread);
        graph.preprocessBedFile(options.bedFile,
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
        graph = std::make_unique<ParallelGraph>(options.numThread);
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
