#include <iostream>
#include <string>
#include "BayesRRm.h"
#include "DenseBayesRRmz.hpp"
#include "data.hpp"
#include "options.hpp"
#include "eigensparsedata.h"
#include "SparseBayesRRG.hpp"
#include "raggedsparsedata.h"
#include "tbb/task_scheduler_init.h"
#include "preprocessgraph.h"

using namespace std;

void processDenseData(Options opt) {
    Data data;

    // Read in the data for every possible option
    data.readFamFile(opt.bedFile + ".fam");
    data.readBimFile(opt.bedFile + ".bim");

    // RAM solution (analysisType = RAMBayes)
    if (opt.analysisType == "RAMBayes" && ( opt.bayesType == "bayes" || opt.bayesType == "bayesMmap" || opt.bayesType == "horseshoe")) {

        clock_t start = clock();

        // Read phenotype file and bed file for the option specified
        data.readPhenotypeFile(opt.phenotypeFile);
        data.readBedFile_noMPI(opt.bedFile+".bed");

        // Option bayesType="bayesMmap" is going to be deprecated
        if (opt.bayesType == "bayesMmap" || opt.bayesType == "bayes"){
            BayesRRm analysis(data, opt, sysconf(_SC_PAGE_SIZE));
            analysis.runGibbs();
        } else if (opt.bayesType == "horseshoe") {
            //TODO Finish horseshoe
        } else if (opt.bayesType == "bayesW") {
            //TODO Add BayesW
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

            PreprocessGraph graph(opt.numThread);
            std::cout << "Preprocessing with " << opt.numThread << " threads ("
                      << (opt.numThreadSpawned > 0 ? std::to_string(opt.numThreadSpawned) : "auto") << " spawned) and "
                      << opt.preprocessChunks << " columns per thread."
                      << endl;
            graph.preprocessBedFile(opt.bedFile + ".bed",
                                    opt.bedFile + ".ppbed",
                                    opt.bedFile + ".ppbedindex",
                                    opt.compress,
                                    &data,
                                    opt.preprocessChunks);
        } else {
            data.preprocessBedFile(opt.bedFile + ".bed",
                    opt.bedFile + ".ppbed",
                    opt.bedFile + ".ppbedindex",
                    opt.compress);
        }

        clock_t end = clock();
        printf("Finished preprocessing the bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
        cout << endl;
    }else if (opt.analysisType == "PPBayes" || opt.analysisType == "PPAsyncBayes") {
        clock_t start = clock();
        data.readPhenotypeFile(opt.phenotypeFile);
        // Run analysis using mapped data files
        if (opt.compress) {
            cout << "Start reading preprocessed bed file: " << opt.bedFile + ".ppbed" << endl;
            clock_t start_bed = clock();
            data.mapCompressedPreprocessBedFile(opt.bedFile + ".ppbed",
                    opt.bedFile + ".ppbedindex");
            clock_t end = clock();
            printf("Finished reading preprocessed bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
            cout << endl;

            std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
            if (opt.numThreadSpawned > 0)
                taskScheduler = std::make_unique<tbb::task_scheduler_init>(opt.numThreadSpawned);

            DenseBayesRRmz analysis(&data, opt);
            analysis.runGibbs();
            data.unmapCompressedPreprocessedBedFile();
        } else {
            cout << "Start reading preprocessed bed file: " << opt.bedFile + ".ppbed" << endl;
            clock_t start_bed = clock();
            data.mapPreprocessBedFile(opt.bedFile + ".ppbed");
            clock_t end = clock();
            printf("Finished reading preprocessed bed file in %.3f sec.\n", double(end - start_bed) / double(CLOCKS_PER_SEC));
            cout << endl;

            BayesRRm analysis(data, opt, sysconf(_SC_PAGE_SIZE));
            analysis.runGibbs();

            data.unmapPreprocessedBedFile();
            end = clock();
            printf("OVERALL read+compute time = %.3f sec.\n", double(end - start) / double(CLOCKS_PER_SEC));
        }
    }else {
        throw(" Error: Wrong analysis type: " + opt.analysisType);
    }
}

void processSparseData(Options options) {
    if (options.analysisType != "PPBayes" && options.analysisType != "PPAsyncBayes") {
        std::cout << "Error: Wrong analysis type: " << options.analysisType << std::endl;
        return;
    }

    using DataPtr = std::unique_ptr<SparseData>;
    DataPtr data;

    if (options.sparseDataType == "eigen") {
        data = DataPtr(new EigenSparseData);
    } else if (options.sparseDataType == "ragged") {
        data = DataPtr(new RaggedSparseData);
    } else {
        std::cout << "Error: Unsupported --sparse-data argument: " << options.sparseDataType << std::endl;
        return;
    }

    // Read in the data for every possible option
    data->readFamFile(options.bedFile + ".fam");
    data->readBimFile(options.bedFile + ".bim");
    data->readPhenotypeFile(options.phenotypeFile);

    // Read the data in sparse format
    data->readBedFileSparse(options.bedFile + ".bed");

    std::unique_ptr<tbb::task_scheduler_init> taskScheduler { nullptr };
    if (options.numThreadSpawned > 0)
        taskScheduler = std::make_unique<tbb::task_scheduler_init>(options.numThreadSpawned);

    SparseBayesRRG analysis(data.get(), options);
    analysis.runGibbs();
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

        if (opt.sparseData)
            processSparseData(opt);
        else
            processDenseData(opt);
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
