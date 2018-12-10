/*
 * BayesRRm.h
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#ifndef SRC_BAYESRRM_H_
#define SRC_BAYESRRM_H_

#include "LinearReg.h"
#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>

class BayesRRm : public LinearReg
{
    int K;
    const double    v0G     = 0.0001;
    const double    s02G    = 0.0001;
    Eigen::VectorXd cva;

    // Component variables
    VectorXd priorPi;   // prior probabilities for each component
    VectorXd pi;        // mixture probabilities
    VectorXd cVa;       // component-specific variance
    VectorXd logL;      // log likelihood of component
    VectorXd muk;       // mean of k-th component marker effect size
    VectorXd denom;     // temporal variable for computing the inflation of the effect variance for a given non-zero componnet
    int m0;             // total num ber of markes in model
    VectorXd v;         // variable storing the component assignment
    VectorXd cVaI;      // inverse of the component variances

    // Mean and residual variables
    double sigmaG;      // genetic variance

    VectorXd Cx;
    VectorXd components;
    std::string header;

    double NM1;
	double km1;
	double betasqn=0;

public:
    BayesRRm(Data &data, Options &opt, const long memPageSize);
    virtual ~BayesRRm();
    void updateBetas(int marker, const VectorXf &Cx);
    void updateHyper();
    void collectSample();


    std::string& getHeader();

private:
    void init(int K, unsigned int markerCount, unsigned int individualCount);
    VectorXd getSnpData(unsigned int marker) const;
    virtual void printDebugInfo();
    virtual void init_header();
};

#endif /* SRC_BAYESRRM_H_ */
