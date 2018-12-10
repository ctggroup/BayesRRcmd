/*
 * Mcmc.cpp
 *
 *  Created on: 6 Dec 2018
 *      Author: admin
 */

#include "Mcmc.h"
#include"options.hpp"
#include "samplewriter.h"
#include <Eigen/Eigen>


Mcmc::Mcmc(Options &opt):
		seed(opt.seed),
		max_iterations(opt.chainLength),
		thinning(opt.thin),
		burn_in(opt.burnin),
		outputFile(opt.mcmcSampleFile){
	    writer.setFileName(outputFile);
	    showDebug=0;
}

Mcmc::~Mcmc() {
	// TODO Auto-generated destructor stub
}
void Mcmc::runGibbs(){
	writer.open(getHeader());
	std::cout << "Running Gibbs sampling" << endl;
	const auto t1 = std::chrono::high_resolution_clock::now();
	for (unsigned int iteration = 0; iteration < max_iterations; iteration++){
		if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
		            std::cout << "iteration: " << iteration << std::endl;
		 sampleGibbsProp();
		 if (iteration >= burn_in && iteration % thinning == 0){
			 writer.write(sample);
		 }
	}
	 const auto t2 = std::chrono::high_resolution_clock::now();
	 const auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	 std::cout << "duration: " << duration << "s" << std::endl;
}

void Mcmc::setDebugEnabled(bool enabled){ showDebug = enabled; }
bool Mcmc::isDebugEnabled(){ return showDebug; }
