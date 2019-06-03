/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef BAYESW_HPP_
#define BAYESW_HPP_

#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>

// Two structures for ARS
struct pars{
	/* Common parameters for the densities */
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on

	VectorXd mixture_classes; // Vector to store mixture component C_k values

	int used_mixture; //Write the index of the mixture we decide to use

	/* Store the current variables */
	double alpha, sigma_b;

	/* Beta_j - specific variables */
	VectorXd X_j;
	/* Mu-specific variables */
	double sigma_mu;
	/* sigma_b-specific variables */
	double alpha_sigma, beta_sigma;

	/*  of sum(X_j*failure) */
	double sum_failure;

	/* Number of events (sum of failure indicators) */
	double d;

	/* Help variable for storing sqrt(2sigma_b)	 */
	double sqrt_2sigmab;

};

struct pars_alpha{
	VectorXd failure_vector;
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on

	/* Alpha-specific variables */
	double alpha_0, kappa_0;  /*  Prior parameters */

	/* Number of events (sum of failure indicators) */
	double d;
};


class BayesW
{
	Data            &data; // data matrices
	Options         &opt;
	const string    bedFile; // bed file
	const long      memPageSize; // size of memory
	const string    outputFile;
	const int       seed;
	const int       max_iterations;
	const int		burn_in;
	const int       thinning;
	const double	alpha_0  = 0.01;
	const double	kappa_0     = 0.01;
	const double    sigma_mu    = 100;
	const double    alpha_sigma  = 1;
	const double    beta_sigma   = 0.0001;
	const string 	quad_points = opt.quad_points;  // Number of Gaussian quadrature points
	const int 		K = opt.S.size()+1;  //number of mixtures + 0 class

	Distributions_boost dist;

	// The ARS variables
	pars used_data;
	pars_alpha used_data_alpha;

	// Component variables
	VectorXd pi_L;        // mixture probabilities
	VectorXd marginal_likelihoods;      // likelihood for each mixture component
	VectorXd v;         // variable storing the component assignment
	VectorXi components; // Indicator vector stating to which mixture SNP belongs to

	// Linear model variables
	VectorXd theta;		 // Fixed effect sizes
	VectorXd beta;       // effect sizes
	VectorXd epsilon;    // variable containing the residuals
	VectorXd vi;		 // adjusted and exponented epsilon

	VectorXd y;
	VectorXd sum_failure;
	VectorXd sum_failure_fix;

	//Sampled variables (not kept in struct)
	double mu;



	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
	double xl, xr ;			  // Initial left and right (pseudo) extremes
	double xinit[4] = {2.5,3,5,10};     // Initial abscissae
	double *p_xinit = xinit;
	VectorXd new_xinit;


public:
	BayesW(Data &data, Options &opt, const long memPageSize);
	virtual ~BayesW();
	int runGibbs_Gauss(); // where we run Gibbs sampling over the parametrised model


private:
	void init(unsigned int markerCount, unsigned int individualCount, unsigned int fixedCount);
	void sampleMu();
	void sampleTheta(int fix_i);
	void sampleBeta(int marker);
	void sampleAlpha();

};


#endif /* BAYESW_HPP_ */
