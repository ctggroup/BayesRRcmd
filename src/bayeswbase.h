/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef BAYESWBASE_H_
#define BAYESWBASE_H_

#include "analysis.h"
#include "common.h"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>
#include <shared_mutex>

struct BayesWKernel;

struct beta_params {
    double alpha = 0;
    double sigma_b = 0;
    double sum_failure = 0;
    double used_mixture = 0;
};

class BayesWBase : public Analysis
{
protected:
    const string    m_bedFile; // bed file
    const long      m_memPageSize; // size of memory
    const string    m_outputFile;
    const int       m_seed;
    const int       m_max_iterations;
    const int		m_burn_in;
    const int       m_thinning;
    const double	m_alpha_0  = 0.01;
    const double	m_kappa_0     = 0.01;
    const double    m_sigma_mu    = 100;
    const double    m_alpha_sigma  = 1;
    const double    m_beta_sigma   = 0.0001;
    const string 	m_quad_points; // Number of Gaussian quadrature points
    const int 		m_K; //number of mixtures + 0 class

    Distributions_boost m_dist;

    VectorXd m_failure_vector;
    double d = 0; // The number of events

	// Component variables
    VectorXd m_mixture_classes; // Vector to store mixture component C_k values
    VectorXd m_pi_L;        // mixture probabilities
    VectorXd m_v;         // variable storing the component assignment
    VectorXi m_components; // Indicator vector stating to which mixture SNP belongs to

	// Linear model variables
    VectorXd m_theta;		 // Fixed effect sizes
    VectorXd m_beta;       // effect sizes
    std::shared_ptr<VectorXd> m_vi = nullptr;		 // adjusted and exponented epsilon

    VectorXd m_y;
    VectorXd m_sum_failure_fix;

    std::shared_ptr<VectorXd> m_epsilon = nullptr; //Vector for residuals
    double m_alpha = 0;
    double m_mu = 0;
    double m_sigma_b = 0;

    std::vector<double> m_randomNumbers;

    mutable std::shared_mutex m_mutex;

public:
    BayesWBase(const Data *data, const Options *opt, const long m_memPageSize);
    virtual ~BayesWBase();

    int runGibbs(AnalysisGraph* analysis) override; // where we run Gibbs sampling over the parametrised model

    void processColumn(const KernelPtr &kernel) override;

    std::unique_ptr<AsyncResult> processColumnAsync(const KernelPtr &kernel) override;
    void doThreadSafeUpdates(const ConstAsyncResultPtr& result) override;
    void updateGlobal(const KernelPtr& kernel, const ConstAsyncResultPtr &result) override;

protected:
	void init(unsigned int markerCount, unsigned int individualCount, unsigned int fixedCount);
	void sampleMu();
    void sampleTheta(int fix_i);
	void sampleAlpha();

    double gauss_hermite_adaptive_integral(int k, double sigma, string n, const BayesWKernel *kernel);

    virtual void prepareForAnalysis();

    virtual int estimateBeta (const BayesWKernel *kernel, const std::shared_ptr<VectorXd> &epsilon, double *xinit, int ninit, double *xl, double *xr, const beta_params params,
                          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                          int nsamp, double *qcent, double *xcent,
                          int ncent, int *neval) = 0;
};


#endif /* BAYESWBASE_H_ */
