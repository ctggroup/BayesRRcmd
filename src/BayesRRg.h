/*
 * BayesRRm.h
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#ifndef SRC_BAYESRRG_H_
#define SRC_BAYESRRG_H_

#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>

class BayesRRg
{
    Data            &data; 			// data matrices
    Options         &opt;
    const string    bedFile; 		// bed file
    const long      memPageSize;	// size of memory
    const string    outputFile;
    const unsigned int seed;
    const unsigned int max_iterations;
    const unsigned int burn_in;
    const unsigned int thinning;
    const double	sigma0  = 0.0001;
    const double	v0E     = 0.0001;
    const double    s02E    = 0.0001;
    const double    v0G     = 0.0001;
    const double    s02G    = 0.0001;
    Distributions_boost dist;
    bool usePreprocessedData;
    bool showDebug;
    double betasqn;

    // Group component variables
    MatrixXd priorPi; 		// prior probabilities for each component
   	MatrixXd pi; 			// mixture probabilities
    VectorXd cVa; 			// component-specific variance
    VectorXd logL; 			// log likelihood of component
    VectorXd muk;			// mean of k-th component marker effect size
    VectorXd denom;			// temporal variable for computing the inflation of the effect variance for a given non-zero component
    int m0;					// total number of markes in model
   	MatrixXd v; 			// variable storing the component assignment
   	VectorXd cVaI;			// inverse of the component variances
    VectorXd components;	// markers M


    // Mean and residual variables
    double mu;          // mean or intercept
    double sigmaG;     	// genetic variance
    double sigmaE;     	// residuals variance
    VectorXd sigmaGG;	// group specific genetic variance

    // Linear model variables
    VectorXd beta;       // effect sizes
    VectorXd y_tilde;    // variable containing the adjusted residuals to exclude the effects of a given marker
    VectorXd epsilon;    // variable containing the residuals

    VectorXd y;
    VectorXd Cx;
    VectorXd betaAcum;

    double epsilonsum;
    double ytildesum;

public:
    BayesRRg(Data &data, Options &opt, const long memPageSize);
    virtual ~BayesRRg();
    int runGibbs(); // where we run Gibbs sampling over the parametrised model

    void setDebugEnabled(bool enabled) { showDebug = enabled; }
    bool isDebugEnabled() const { return showDebug; }

private:
    void init(int K, unsigned int markerCount, unsigned int individualCount);
    VectorXd getSnpData(unsigned int marker) const;
    void printDebugInfo() const;
};

#endif /* SRC_BAYESRRG_H_ */
