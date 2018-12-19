/*
 * LinearReg.cpp
 *
 *  Created on: 7 Dec 2018
 *      Author: admin
 */

#include "LinearReg.h"
#include "data.hpp"
#include "distributions_boost.hpp"
#include "options.hpp"
#include "parallelalgo.h"

#include <random>

LinearReg::LinearReg(Options &opt,Data &data):
		     Mcmc(opt)
		    , data(data)
    	    , M(data.numIncdSnps)
    	    , N(data.numKeptInds)
    	    , opt(opt)
    	    , bedFile(opt.bedFile + ".bed")
    	    , memPageSize(memPageSize)
    	    , dist(opt.seed)
    	    , usePreprocessedData(opt.analysisType == "PPBayes")
            {
	 // Linear model variables
		      beta = VectorXd(M);           // effect sizes
		      y_tilde = VectorXd(N);    // variable containing the adjusted residuals to exclude the effects of a given marker
		      epsilon = VectorXd(N);    // variable containing the residuals

		      y = VectorXd();
		      y_tilde.setZero();
		      beta.setZero();
		      y = (data.y.cast<double>().array() - data.y.cast<double>().mean());
		      y /= sqrt(y.squaredNorm() / (double(N - 1)));

		     epsilon = (y).array() - mu;
		     sigmaE = epsilon.squaredNorm() / M * 0.5;
		     markerI.resize(M);
		     std::iota(markerI.begin(), markerI.end(), 0);
		     mu=0;


}

LinearReg::~LinearReg() {
	// TODO Auto-generated destructor stub


}

void LinearReg::sampleGibbsProp(){

	const double sigmaEpsilon = parallelStepAndSumEpsilon(epsilon, mu);
	parallelStepMuEpsilon(mu,epsilon, sigmaEpsilon, double(N), sigmaE, dist);
    std::random_shuffle(markerI.begin(), markerI.end());
    for (unsigned int j = 0; j < M; j++) {


    	const auto marker = markerI[j];
    	const VectorXf &Cx = data.mappedZ.col(marker);
        // Now y_tilde = Y-mu - X * beta + X.col(marker) * beta(marker)_old
    	double beta_old=beta(marker);
    	if(beta_old==0) //if beta=0 then there is no need to update y_tilde
    	 parallelUpdateYTilde(y_tilde, epsilon, Cx, beta(marker));
    	else
    		y_tilde=epsilon;
    	updateBetas(marker,Cx);
    	 // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
    	if(beta(marker)!=0){
    		parallelUpdateEpsilon(epsilon, y_tilde, Cx, beta(marker));
    		betasqn+=pow(beta(marker),2)-pow(beta_old,2);
    	} //no need to update epsilon if beta==0

    }
    updateHyper();

    const double epsilonSqNorm = parallelSquaredNorm(epsilon);
    sigmaE = dist.inv_scaled_chisq_rng(v0E + N, (epsilonSqNorm + v0E * s02E) / (v0E + N));
    std::cout<<"sigmaE "<<sigmaE<<"\n";

    if (showDebug)
        printDebugInfo();
    collectSample();
}



