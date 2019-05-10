#ifndef options_hpp
#define options_hpp

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <string>
#include <limits.h>
#include <boost/format.hpp>
#include "gadgets.hpp"
#include <Eigen/Eigen>
#include "common.h"

using namespace std;
using namespace boost;

const unsigned Megabase = 1e6;

class Options {
public:
    unsigned chainLength;
    unsigned burnin;
    unsigned seed;
    unsigned numThread = 0; // Default to tbb::flow::unlimited
    int numThreadSpawned = 0; // Default to 0, let TBB do its thing
    size_t decompressionNodeConcurrency = 0;
    size_t decompressionTokens = 40;
    size_t analysisNodeConcurrency = 0;
    size_t analysisTokens = 20;
    unsigned preprocessChunks = 1;
    unsigned thin;  // save every this th sampled value in MCMC
    vector<float> S;    //variance components

    unsigned int numGroups;
    Eigen::MatrixXd mS;
    string groupFile;
   
    string title;
    string analysisType;
    string bayesType;
    string phenotypeFile;
    string bedFile;
    string mcmcSampleFile;
    string optionFile;
    bool compress = false;
    DataType dataType = DataType::Dense;
    string iterLogFile;
    bool iterLog = false;
    string colLogFile;
    bool colLog =false;

    Options(){
        chainLength             = 10000;
        burnin                  = 5000;
        seed                    = static_cast<unsigned int>(std::time(0));
        numThread               = 0;
        numThreadSpawned        = 0;
        decompressionNodeConcurrency = 0;
        decompressionTokens     = 40;
        analysisNodeConcurrency = 0;
        analysisTokens          = 20;
        preprocessChunks        = 1;
        thin                    = 5;

        S.resize(3);
        S[0]                    = 0.01;
        S[1]                    = 0.001;
        S[2]                    = 0.0001;

        title                   = "brr";
        analysisType            = "Bayes";
        bayesType               = "C";
        phenotypeFile           = "";
        bedFile                 = "";
        mcmcSampleFile          = "bayesOutput.csv";
        optionFile				= "";
        numGroups				=2;
        dataType                = DataType::Dense;
    }

    void inputOptions(const int argc, const char* argv[]);

private:
    void readFile(const string &file);
    void makeTitle(void);
    void seedEngine(void);
};

#endif /* options_hpp */
