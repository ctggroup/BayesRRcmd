/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef BAYESWBASE_H_
#define BAYESWBASE_H_

#include "bayeswcommon.h"
#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>

#if 0
// Two structures for ARS
struct pars{
	/* Common parameters for the densities */
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on

	VectorXd mixture_classes; // Vector to store mixture component C_k values

	int used_mixture; //Write the index of the mixture we decide to use

	/* Store the current variables */
	double alpha;

	/* Beta_j - specific variables */
	VectorXd X_j;

	/*  of sum(X_j*failure) */
	double sum_failure;

	/* Mu-specific variables */
	double sigma_mu;
	/* sigma_b-specific variables */
	double alpha_sigma, beta_sigma;

	/* Number of events (sum of failure indicators) */
	double d;

	/* Help variable for storing sqrt(2sigma_b)	 */
	double sqrt_2sigmab;

};

struct pars_beta_sparse{
	/* Common parameters for the densities */

	VectorXd mixture_classes; // Vector to store mixture component C_k values
    // assigned in init

	int used_mixture; //Write the index of the mixture we decide to use
    // assigned in sampleBeta, then used within the same scope

	/* Store the current variables */
	double alpha, sigma_b;
    // alpha - initialised in init, updated in sampleAlpha - in the markov chain, after sampleBeta
    // sigma_b - initialised in init, updated in the markov chain after sampleBeta

	/* Beta_j - specific variables */
	double vi_0, vi_1, vi_2; // Sums of vi elements
    // assigned in sampleBeta, then used within the same scope

	// Mean, std dev and their ratio for snp j
	double mean, sd, mean_sd_ratio;
    // mean - assined in sampleBeta but not used
    // sd - assigned in sampleBeta, used in sampleBeta and beta_dens
    // mean_sd_ratio - assigned in sampleBeta, used in sampleBeta and beta_dens


	/*  of sum(X_j*failure) */
	double sum_failure;
    // assigned in sampleBeta, used in gauss_hermite_adaptive_integral and alpha_dens
};

struct pars_alpha{
	VectorXd failure_vector;
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on

	/* Alpha-specific variables */
	double alpha_0, kappa_0;  /*  Prior parameters */

	/* Number of events (sum of failure indicators) */
	double d;
};

#endif

struct gh_params {
    virtual double exponent_sum() const = 0;
    virtual double integrand_adaptive(double s,double alpha, double dj, double sqrt_2Ck_sigmab) const = 0;
};

class BayesWBase
{
protected:
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
	struct pars used_data;
	struct pars_beta_sparse used_data_beta;
	struct pars_alpha used_data_alpha;

	// Component variables
	VectorXd pi_L;        // mixture probabilities
	VectorXd marginal_likelihoods;      // likelihood for each mixture component
	VectorXd v;         // variable storing the component assignment
	VectorXi components; // Indicator vector stating to which mixture SNP belongs to

	// Linear model variables
	VectorXd theta;		 // Fixed effect sizes
	VectorXd beta;       // effect sizes
	VectorXd vi;		 // adjusted and exponented epsilon

	VectorXd y;
	VectorXd sum_failure;
	VectorXd sum_failure_fix;

	VectorXd epsilon; //Vector for residuals
	//Sampled variables (not kept in struct)
	double mu;



public:
        BayesWBase(Data &data, Options &opt, const long memPageSize);
        virtual ~BayesWBase();
	int runGibbs_Gauss(); // where we run Gibbs sampling over the parametrised model


protected:
	void init(unsigned int markerCount, unsigned int individualCount, unsigned int fixedCount);
	void sampleMu();
	void sampleTheta(int fix_i);
    virtual void sampleBeta(int marker);
	void sampleAlpha();

    virtual double calculateSumFailure(int marker) = 0;

    void marginal_likelihood_vec_calc(VectorXd prior_prob, VectorXd &post_marginals, string n, const gh_params *params);
    double gauss_hermite_adaptive_integral(int k, double sigma, string n, const gh_params *params);

    virtual void preEstimateResidualUpdate(int marker) = 0;

    virtual std::unique_ptr<gh_params> gaussHermiteParameters(int marker) = 0;

    virtual int estimateBeta (double *xinit, int ninit, double *xl, double *xr,
                          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                          int nsamp, double *qcent, double *xcent,
                          int ncent, int *neval) = 0;

    virtual void postEstimateResidualUpdate(int marker) = 0;
};


#endif /* BAYESWBASE_H_ */
