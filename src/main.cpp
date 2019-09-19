#include <iostream>

#include "analysisrunner.h"
#include "gadgets.hpp"
#include "options.hpp"
#include "version.h"

int main(int argc, const char * argv[]) {

    using namespace std;

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
        cout<<"Software version info: " << GIT_COMMIT << endl;
	cout<<"Seed used: " << options.seed << endl;
        if (!AnalysisRunner::run(options))
            exit(1);
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
