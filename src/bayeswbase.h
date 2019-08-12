/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef BAYESWBASE_H_
#define BAYESWBASE_H_

#include "analysis.h"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>

struct BayesWKernel;
struct Kernel;
struct Marker;
struct MarkerBuilder;

struct beta_params {
    double alpha = 0;
    double sigma_b = 0;
    double sum_failure = 0;
    double used_mixture = 0;
};

class BayesWBase : public Analysis
{
protected:
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
    const string 	quad_points; // Number of Gaussian quadrature points
    const int 		K; //number of mixtures + 0 class

	Distributions_boost dist;

    VectorXd failure_vector;
    double d = 0; // The number of events

	// Component variables
    VectorXd mixture_classes; // Vector to store mixture component C_k values
	VectorXd pi_L;        // mixture probabilities
	VectorXd marginal_likelihoods;      // likelihood for each mixture component
	VectorXd v;         // variable storing the component assignment
	VectorXi components; // Indicator vector stating to which mixture SNP belongs to

	// Linear model variables
	VectorXd theta;		 // Fixed effect sizes
	VectorXd beta;       // effect sizes
	VectorXd vi;		 // adjusted and exponented epsilon

    VectorXd y;
	VectorXd sum_failure_fix;

    VectorXd epsilon; //Vector for residuals
    double alpha = 0;
    double mu = 0;
    double sigma_b = 0;



public:
    BayesWBase(const Data *data, const Options *opt, const long memPageSize);
    virtual ~BayesWBase();

    int runGibbs(AnalysisGraph* analysis) override; // where we run Gibbs sampling over the parametrised model

    void processColumn(Kernel *kernel) override;
    std::unique_ptr<AsyncResult> processColumnAsync(Kernel *kernel) override;

    void updateGlobal(Kernel *kernel, const double beta_old, const double beta, const VectorXd &deltaEps) override;

protected:
	void init(unsigned int markerCount, unsigned int individualCount, unsigned int fixedCount);
	void sampleMu();
    void sampleTheta(int fix_i);
	void sampleAlpha();

    double gauss_hermite_adaptive_integral(int k, double sigma, string n, const BayesWKernel *kernel);

    virtual int estimateBeta (const BayesWKernel *kernel, double *xinit, int ninit, double *xl, double *xr, const beta_params params,
                          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                          int nsamp, double *qcent, double *xcent,
                          int ncent, int *neval) = 0;
};


#endif /* BAYESWBASE_H_ */
