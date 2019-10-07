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
    static MatrixXd parseVarianceComponents(const std::string &arg);

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
    Eigen::MatrixXd S;    //variance components

    unsigned int numGroups;
    string groupFile;
    string failureFile;
    string bayesW_version;
    string quad_points;
    string fixedFile;
    unsigned int fixedEffectNumber;

    string title;
    AnalysisType analysisType = AnalysisType::Unknown;
    string phenotypeFile;
    string dataFile;
    InputType inputType = InputType::Unknown;
    string mcmcSampleFile;
    string optionFile;
    bool compress = false;
    PreprocessDataType preprocessDataType = PreprocessDataType::Dense;
    string iterLogFile;
    bool iterLog = false;
    string colLogFile;
    bool colLog =false;
    bool useMarkerCache = false;


    double v0E  = 0.0001;
    double s02E = 0.0001;
    double v0G  = 0.0001;
    double s02G = 0.0001;
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

        S.resize(1, 3);
        S(0,0)                  = 0.01;
        S(0,1)                  = 0.001;
        S(0,2)                  = 0.0001;

        title                   = "brr";
        analysisType            = AnalysisType::Unknown;
        dataFile                = "";
        inputType               = InputType::Unknown;
        fixedFile 				= "";
        phenotypeFile           = "";
        mcmcSampleFile          = "bayesOutput.csv";
        optionFile				= "";
        numGroups				=2;
        preprocessDataType      = PreprocessDataType::Dense;

        bayesW_version		= "marginal";
        fixedEffectNumber       = 0;
    }

    void inputOptions(const int argc, const char* argv[]);

private:
    void readFile(const string &file);
    void makeTitle(void);
    void seedEngine(void);
};

#endif /* options_hpp */
