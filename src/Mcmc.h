/*
 * Mcmc.h
 *
 *  Created on: 6 Dec 2018
 *      Author: Daniel Trejo Banos
 */

#ifndef MCMC_H_
#define MCMC_H_
#include"options.hpp"
#include "samplewriter.h"
#include <Eigen/Eigen>

class Mcmc {
	public:

	const string  outputFile;
    const int     seed;
	const int 	  max_iterations;
    const int	  burn_in;
	const int 	  thinning;
	SampleWriter  writer;
	VectorXd	  sample;
    bool showDebug;

	Mcmc(Options &opt);
	virtual ~Mcmc();
	void runGibbs();
	virtual void sampleGibbsProp()=0;
	virtual std::string& getHeader()=0;
	virtual  void printDebugInfo()=0;
	virtual  void collectSample()=0;
	void setDebugEnabled(bool enabled);
    bool isDebugEnabled();
};

#endif /* MCMC_H_ */
