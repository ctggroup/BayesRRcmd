#include <iostream>
#include <string>
#include "BayesRRm.h"
#include "data.hpp"
#include "options.hpp"



using namespace std;

int main(int argc, const char * argv[]) {

    if (argc < 2){
        cerr << " \nDid you forget to give the input parameters?\n" << endl;
        exit(1);
    }

    try {
        Options opt;
        opt.inputOptions(argc, argv);

        Data data;

        // Read in the data for every possible option
        data.readFamFile(opt.bedFile + ".fam");
        data.readBimFile(opt.bedFile + ".bim");



        if (opt.bedToSparse == true) {
            
            data.readPhenotypeFile(opt.phenotypeFile);
            BayesRRm analysis(data, opt, sysconf(_SC_PAGE_SIZE));
            analysis.write_sparse_data_files();

        } else if (opt.bayesType == "bayesMPI" && opt.analysisType == "RAM") {
            // Read phenotype file
            data.readPhenotypeFile(opt.phenotypeFile);
            //EO to remove
            //data.readBedFile_noMPI(opt.bedFile+".bed");
            BayesRRm analysis(data, opt, sysconf(_SC_PAGE_SIZE));

            if (opt.markerBlocksFile != "") {
                data.readMarkerBlocksFile(opt.markerBlocksFile);
            } else {
                cout << "No request to read a block markers definition file." << endl;
            }
            analysis.runMpiGibbs();

        }
        // Pre-processing the data (centering and scaling)
        
        else {
            throw(" Error: Wrong analysis requested: " + opt.analysisType + " + " + opt.bayesType);
        }

        //#endif

    }
        
    catch (const string &err_msg) {
        cerr << "\n" << err_msg << endl;
    }
    catch (const char *err_msg) {
        cerr << "\n" << err_msg << endl;
    }

    return 0;
}
