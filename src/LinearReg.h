/*
 * LinearReg.h
 *
 *  Created on: 7 Dec 2018
 *      Author: admin
 */

#ifndef LINEARREG_H_
#define LINEARREG_H_

#include "Mcmc.h"
#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

class LinearReg : public Mcmc {


public:
	    Distributions_boost dist;
	    bool usePreprocessedData;
		const int M;
		const int N;
	    Data            &data; // data matrices
	    Options         &opt;
	    const string    bedFile; // bed file
	    const long      memPageSize; // size of memory
	    const double	sigma0  = 0.0001;
	    const double	v0E     = 0.0001;
	    const double    s02E    = 0.0001;
	    // Linear model variables
	     VectorXd beta;       // effect sizes
	     VectorXd y_tilde;    // variable containing the adjusted residuals to exclude the effects of a given marker
	     VectorXd epsilon;    // variable containing the residuals

	     VectorXd y;
	     double sigmaE;      // residuals variance
	     double mu;          // mean or intercept
	     std::vector<unsigned long int> markerI;
	 	double betasqn=0;


	LinearReg(Options &opt,Data &data);
	virtual ~LinearReg();
    void sampleGibbsProp();
    virtual void updateBetas(int marker, const VectorXf &Cx)=0;
    virtual void updateHyper()=0;

};

#endif /* LINEARREG_H_ */
