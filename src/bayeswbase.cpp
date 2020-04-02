/*
 * BayesW.cpp
 *
 *  Created on: 26 Nov 2018
 *  Author: Sven Erik Ojavee
 *  Last changes: 22 Feb 2019
 */

#include "analysisgraph.hpp"
#include "data.hpp"
#include "distributions_boost.hpp"
//#include "concurrentqueue.h"
#include "options.hpp"
#include "bayeswbase.h"
#include "bayeswkernel.h"
#include "BayesW_arms.h"
#include "markerbuilder.h"
#include "samplewriter.h"

#include <chrono>
#include <numeric>
#include <random>

/* Pre-calculate used constants */
#define PI 3.14159
#define PI2 6.283185
#define sqrtPI 1.77245385090552
#define EuMasc 0.577215664901532

BayesWBase::BayesWBase(const Data *data, const Options *opt, const long memPageSize)
: Analysis(data, opt)
, m_seed(opt->seed)
, m_memPageSize(memPageSize)
, m_max_iterations(opt->chainLength)
, m_thinning(opt->thin)
, m_burn_in(opt->burnin)
, m_outputFile(opt->mcmcSampleFile)
, m_bedFile(opt->dataFile)
, m_dist(opt->seed)
, m_quad_points(opt->quad_points)
, m_K(opt->S.size() + 1)
{

}

BayesWBase::~BayesWBase() = default;

namespace  {

/* Function to check if ARS resulted with error*/
inline void errorCheck(int err){
	if(err>0){
		cout << "Error code = " << err << endl;
		exit(1);
	}
}

struct mu_params {
    double alpha = 0;
    double d = 0;
    VectorXd epsilon;
    double sigma_mu = 0;
};

/* Function for the log density of mu */
inline double mu_dens(double x, void *norm_data)
/* We are sampling mu (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
    mu_params p = *(static_cast<mu_params *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
    y = - p.alpha * x * p.d - (( p.epsilon.array()  - x) * p.alpha - EuMasc).exp().sum() - x*x/(2*p.sigma_mu);
	return y;
};

struct theta_params {
    double alpha = 0;
    double sum_failure = 0;
    VectorXd epsilon;
    VectorXd X_j;
    double sigma_mu = 0;
};

/* Function for the log density of some "fixed" covariate effect */
inline double theta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;
	/* In C++ we need to do a static cast for the void data */
    theta_params p = *(static_cast<theta_params *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = - p.alpha * x * p.sum_failure - (((p.epsilon -  p.X_j * x)* p.alpha).array() - EuMasc).exp().sum() - x*x/(2*p.sigma_mu); // Prior is the same currently for intercepts and fixed effects
	return y;
};

struct alpha_params {
    double alpha_0 = 0;
    double d = 0;
    VectorXd epsilon;
    VectorXd failure_vector;
    double kappa_0 = 0;
};

/* Function for the log density of alpha */
inline double alpha_dens(double x, void *norm_data)
/* We are sampling alpha (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
    alpha_params p = *(static_cast<alpha_params *>(norm_data));
	y = (p.alpha_0 + p.d - 1) * log(x) + x * ((p.epsilon.array() * p.failure_vector.array()).sum() - p.kappa_0) -
			((p.epsilon * x).array() - EuMasc).exp().sum() ;
	return y;
};

}


//Calculate the value of the integral using Adaptive Gauss-Hermite quadrature
//Let's assume that mu is always 0 for speed
double BayesWBase::gauss_hermite_adaptive_integral(int k, double sigma, string n, const BayesWKernel *kernel){
    assert(kernel);

    double temp = 0;
    double sqrt_2ck_sigma = sqrt(2*m_mixture_classes(k)*m_sigma_b);

    if(n == "3"){
        double x1,x2;
        double w1,w2,w3;

        x1 = 1.2247448713916;
        x2 = -x1;

        w1 = 1.3239311752136;
        w2 = w1;

        w3 = 1.1816359006037;

        x1 = sigma*x1;
        x2 = sigma*x2;

        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3;
    }
    // n=5
    else if(n == "5"){
        double x1,x2,x3,x4;//x5;
        double w1,w2,w3,w4,w5; //These are adjusted weights

        x1 = 2.0201828704561;
        x2 = -x1;
        w1 = 1.181488625536;
        w2 = w1;

        x3 = 0.95857246461382;
        x4 = -x3;
        w3 = 0.98658099675143;
        w4 = w3;

        //	x5 = 0.0;
        w5 = 0.94530872048294;

        x1 = sigma*x1;
        x2 = sigma*x2;
        x3 = sigma*x3;
        x4 = sigma*x4;
        //x5 = sigma*x5;

        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3 * kernel->integrand_adaptive(x3,m_alpha,sqrt_2ck_sigma) +
                w4 * kernel->integrand_adaptive(x4,m_alpha,sqrt_2ck_sigma) +
                w5;
    }else if(n == "7"){
        double x1,x2,x3,x4,x5,x6;
        double w1,w2,w3,w4,w5,w6,w7; //These are adjusted weights

        x1 = 2.6519613568352;
        x2 = -x1;
        w1 = 1.1013307296103;
        w2 = w1;

        x3 = 1.6735516287675;
        x4 = -x3;
        w3 = 0.8971846002252;
        w4 = w3;

        x5 = 0.81628788285897;
        x6 = -x5;
        w5 = 0.8286873032836;
        w6 = w5;

        w7 = 0.81026461755681;

        x1 = sigma*x1;
        x2 = sigma*x2;
        x3 = sigma*x3;
        x4 = sigma*x4;
        x5 = sigma*x5;
        x6 = sigma*x6;

        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3 * kernel->integrand_adaptive(x3,m_alpha,sqrt_2ck_sigma) +
                w4 * kernel->integrand_adaptive(x4,m_alpha,sqrt_2ck_sigma) +
                w5 * kernel->integrand_adaptive(x5,m_alpha,sqrt_2ck_sigma) +
                w6 * kernel->integrand_adaptive(x6,m_alpha,sqrt_2ck_sigma) +
                w7;
    }else if(n == "9"){
        double x1,x2,x3,x4,x5,x6,x7,x8;//x9;
        double w1,w2,w3,w4,w5,w6,w7,w8,w9; //These are adjusted weights

        x1 = 3.1909932017815;
        x2 = -x1;
        w1 = 1.0470035809767;
        w2 = w1;

        x3 = 2.2665805845318;
        x4 = -x3;
        w3 = 0.84175270147867;
        w4 = w3;

        x5 = 1.4685532892167;
        x6 = -x3;
        w5 = 0.7646081250946;
        w6 = w5;

        x7 = 0.72355101875284;
        x8 = -x7;
        w7 = 0.73030245274509;
        w8 = w7;

        //x9 = 0.0;
        w9 = 0.72023521560605;

        x1 = sigma*x1;
        x2 = sigma*x2;
        x3 = sigma*x3;
        x4 = sigma*x4;
        x5 = sigma*x5;
        x6 = sigma*x6;
        x7 = sigma*x7;
        x8 = sigma*x8;


        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3 * kernel->integrand_adaptive(x3,m_alpha,sqrt_2ck_sigma) +
                w4 * kernel->integrand_adaptive(x4,m_alpha,sqrt_2ck_sigma) +
                w5 * kernel->integrand_adaptive(x5,m_alpha,sqrt_2ck_sigma) +
                w6 * kernel->integrand_adaptive(x6,m_alpha,sqrt_2ck_sigma) +
                w7 * kernel->integrand_adaptive(x7,m_alpha,sqrt_2ck_sigma) +
                w8 * kernel->integrand_adaptive(x8,m_alpha,sqrt_2ck_sigma) +
                w9;
    }else if(n == "11"){
        double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;//,x11;
        double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11; //These are adjusted weights

        x1 = 3.6684708465596;
        x2 = -x1;
        w1 = 1.0065267861724;
        w2 = w1;

        x3 = 2.7832900997817;
        x4 = -x3;
        w3 = 0.802516868851;
        w4 = w3;

        x5 = 2.0259480158258;
        x6 = -x3;
        w5 = 0.721953624728;
        w6 = w5;

        x7 = 1.3265570844949;
        x8 = -x7;
        w7 = 0.6812118810667;
        w8 = w7;

        x9 = 0.6568095668821;
        x10 = -x9;
        w9 = 0.66096041944096;
        w10 = w9;

        //x11 = 0.0;
        w11 = 0.65475928691459;

        x1 = sigma*x1;
        x2 = sigma*x2;
        x3 = sigma*x3;
        x4 = sigma*x4;
        x5 = sigma*x5;
        x6 = sigma*x6;
        x7 = sigma*x7;
        x8 = sigma*x8;
        x9 = sigma*x9;
        x10 = sigma*x10;
        //	x11 = sigma*x11;

        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3 * kernel->integrand_adaptive(x3,m_alpha,sqrt_2ck_sigma) +
                w4 * kernel->integrand_adaptive(x4,m_alpha,sqrt_2ck_sigma) +
                w5 * kernel->integrand_adaptive(x5,m_alpha,sqrt_2ck_sigma) +
                w6 * kernel->integrand_adaptive(x6,m_alpha,sqrt_2ck_sigma) +
                w7 * kernel->integrand_adaptive(x7,m_alpha,sqrt_2ck_sigma) +
                w8 * kernel->integrand_adaptive(x8,m_alpha,sqrt_2ck_sigma) +
                w9 * kernel->integrand_adaptive(x9,m_alpha,sqrt_2ck_sigma) +
                w10 * kernel->integrand_adaptive(x10,m_alpha,sqrt_2ck_sigma) +
                w11;
    }else if(n == "13"){
        double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12;
        double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13; //These are adjusted weights

        x1 = 4.1013375961786;
        x2 = -x1;
        w1 = 0.97458039564;
        w2 = w1;

        x3 = 3.2466089783724;
        x4 = -x3;
        w3 = 0.7725808233517;
        w4 = w3;

        x5 = 2.5197356856782;
        x6 = -x3;
        w5 = 0.6906180348378;
        w6 = w5;

        x7 = 1.8531076516015;
        x8 = -x7;
        w7 = 0.6467594633158;
        w8 = w7;

        x9 = 1.2200550365908;
        x10 = -x9;
        w9 = 0.6217160552868;
        w10 = w9;

        x11 = 0.60576387917106;
        x12 = -x11;
        w11 = 0.60852958370332;
        w12 = w11;

        //x13 = 0.0;
        w13 = 0.60439318792116;

        x1 = sigma*x1;
        x2 = sigma*x2;
        x3 = sigma*x3;
        x4 = sigma*x4;
        x5 = sigma*x5;
        x6 = sigma*x6;
        x7 = sigma*x7;
        x8 = sigma*x8;
        x9 = sigma*x9;
        x10 = sigma*x10;
        x11 = sigma*x11;
        x12 = sigma*x12;

        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3 * kernel->integrand_adaptive(x3,m_alpha,sqrt_2ck_sigma) +
                w4 * kernel->integrand_adaptive(x4,m_alpha,sqrt_2ck_sigma) +
                w5 * kernel->integrand_adaptive(x5,m_alpha,sqrt_2ck_sigma) +
                w6 * kernel->integrand_adaptive(x6,m_alpha,sqrt_2ck_sigma) +
                w7 * kernel->integrand_adaptive(x7,m_alpha,sqrt_2ck_sigma) +
                w8 * kernel->integrand_adaptive(x8,m_alpha,sqrt_2ck_sigma) +
                w9 * kernel->integrand_adaptive(x9,m_alpha,sqrt_2ck_sigma) +
                w10 * kernel->integrand_adaptive(x10,m_alpha,sqrt_2ck_sigma) +
                w11 * kernel->integrand_adaptive(x11,m_alpha,sqrt_2ck_sigma) +
                w12 * kernel->integrand_adaptive(x12,m_alpha,sqrt_2ck_sigma) +
                w13;
    }else if(n == "15"){
        double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14;//,x11;
        double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15; //These are adjusted weights

        x1 = 4.4999907073094;
        x2 = -x1;
        w1 = 0.94836897082761;
        w2 = w1;

        x3 = 3.6699503734045;
        x4 = -x3;
        w3 = 0.7486073660169;
        w4 = w3;

        x5 = 2.9671669279056;
        x6 = -x3;
        w5 = 0.666166005109;
        w6 = w5;

        x7 = 2.3257324861739;
        x8 = -x7;
        w7 = 0.620662603527;
        w8 = w7;

        x9 = 1.7199925751865;
        x10 = -x9;
        w9 = 0.5930274497642;
        w10 = w9;

        x11 = 1.1361155852109;
        x12 = -x11;
        w11 = 0.5761933502835;
        w12 = w11;

        x13 = 0.5650695832556;
        x14 = -x13;
        w13 = 0.5670211534466;
        w14 = w13;

        //x15 = 0.0;
        w15 = 0.56410030872642;

        x1 = sigma*x1;
        x2 = sigma*x2;
        x3 = sigma*x3;
        x4 = sigma*x4;
        x5 = sigma*x5;
        x6 = sigma*x6;
        x7 = sigma*x7;
        x8 = sigma*x8;
        x9 = sigma*x9;
        x10 = sigma*x10;
        x11 = sigma*x11;
        x12 = sigma*x12;
        x13 = sigma*x13;
        x14 = sigma*x14;

        temp = w1 * kernel->integrand_adaptive(x1,m_alpha,sqrt_2ck_sigma) +
                w2 * kernel->integrand_adaptive(x2,m_alpha,sqrt_2ck_sigma) +
                w3 * kernel->integrand_adaptive(x3,m_alpha,sqrt_2ck_sigma) +
                w4 * kernel->integrand_adaptive(x4,m_alpha,sqrt_2ck_sigma) +
                w5 * kernel->integrand_adaptive(x5,m_alpha,sqrt_2ck_sigma) +
                w6 * kernel->integrand_adaptive(x6,m_alpha,sqrt_2ck_sigma) +
                w7 * kernel->integrand_adaptive(x7,m_alpha,sqrt_2ck_sigma) +
                w8 * kernel->integrand_adaptive(x8,m_alpha,sqrt_2ck_sigma) +
                w9 * kernel->integrand_adaptive(x9,m_alpha,sqrt_2ck_sigma) +
                w10 * kernel->integrand_adaptive(x10,m_alpha,sqrt_2ck_sigma) +
                w11 * kernel->integrand_adaptive(x11,m_alpha,sqrt_2ck_sigma) +
                w12 * kernel->integrand_adaptive(x12,m_alpha,sqrt_2ck_sigma) +
                w13 * kernel->integrand_adaptive(x13,m_alpha,sqrt_2ck_sigma) +
                w14 * kernel->integrand_adaptive(x14,m_alpha,sqrt_2ck_sigma) +
                w15;
    }else{
        cout << "Possible number of quad_points = 3,5,7,9,11,13,15" << endl;
        exit(1);
    }

    return sigma*temp;
}

void BayesWBase::prepareForAnalysis()
{
    // Generate the random numbers required for this iteration. The random
    // number engine is not thread safe, so generate them up front to avoid
    // having to use a mutex.
    std::generate(m_randomNumbers.begin(), m_randomNumbers.end(), [&dist = m_dist]() {
        return dist.unif_rng();
    });
}

void BayesWBase::init(unsigned int markerCount, unsigned int individualCount, unsigned int fixedCount)
{
	// Component variables
    m_pi_L = VectorXd(m_K);           		 // prior mixture probabilities
    m_v = VectorXd(m_K);            		 // vector storing the component assignment

	// Linear model variables
    m_beta = VectorXd(markerCount);           // effect sizes
    m_theta = VectorXd(fixedCount);

    m_sum_failure_fix = VectorXd(fixedCount); // Vector to sum fixed vector * failure vector per fixed effect

	//phenotype vector
    m_y = VectorXd();

    m_vi = std::make_shared<VectorXd>(individualCount);

	// Init the working variables
    const int km1 = m_K - 1;

	//vector with component class for each marker
    m_components = VectorXi(markerCount);
    m_components.setZero();

	//set priors for pi parameters
	//Give all mixtures (except 0 class) equal initial probabilities
    m_pi_L(0) = 0.99;
    m_pi_L.segment(1,km1).setConstant((1-m_pi_L(0))/km1);

    m_beta.setZero();
    m_theta.setZero();

	//initialize epsilon vector as the phenotype vector
    m_y = m_data->y.cast<double>().array();
    m_mu = m_y.mean();       // mean or intercept

	// Initialize the variables in structures
	//Save variance classes
    m_mixture_classes.resize(km1); //The future solution


	for(int i=0;i<(km1);i++){
        m_mixture_classes(i) = m_opt->S.row(0)[i];   //Save the mixture data (C_k)
	}

    //Store the vector of failures
    m_failure_vector = m_data->fail.cast<double>();
    // Save the number of events
    d = m_failure_vector.array().sum();

    double denominator = (6 * ((m_y.array() - m_mu).square()).sum()/(m_y.size()-1));
    m_alpha = PI/sqrt(denominator);    // The shape parameter initial value

    m_epsilon = std::make_shared<VectorXd>(m_y.array() - m_mu); // Initially, all the BETA elements are set to 0, XBeta = 0

    m_sigma_b = PI2/ (6 * pow(m_alpha,2) * markerCount ) ;

	// Save the sum(X_j*failure) for each j
	//Previous preprocessed version for reading columns
	//for(int marker=0; marker<M; marker++){
    //	sum_failure(marker) = ((m_data->mappedZ.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	//}

	//If there are fixed effects, find the same values for them
	if(fixedCount > 0){
		for(int fix_i=0; fix_i < fixedCount; fix_i++){
            m_sum_failure_fix(fix_i) = ((m_data->X.col(fix_i).cast<double>()).array() * m_failure_vector.array()).sum();
		}
	}

    m_randomNumbers.resize(markerCount);
}
// Function for sampling intercept (mu)
void BayesWBase::sampleMu(){
	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
    double xinit[4] = {0.995*m_mu, m_mu,  1.005*m_mu, 1.01*m_mu};     // Initial abscissae
	double *p_xinit = xinit;

	double xl = 2;
	double xr = 5;   //xl and xr and the maximum and minimum values between which we sample

    mu_params params;
    params.alpha = m_alpha;
    params.d = d;
    params.epsilon = m_epsilon->array() + m_mu; // we add to epsilon =Y+mu-X*beta
    params.sigma_mu = m_sigma_mu;

	// Use ARS to sample mu (with density mu_dens, using parameters from used_data)
    err = arms(xinit,ninit,&xl,&xr,mu_dens,&params,&convex,
			npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

	errorCheck(err); // If there is error, stop the program
    m_mu = xsamp[0];   // Save the sampled value
    *m_epsilon = params.epsilon.array() - m_mu;// we substract again now epsilon =Y-mu-X*beta
}

// Function for sampling fixed effect (theta_i)
void BayesWBase::sampleTheta(int fix_i){
	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
    double xinit[4] = {m_theta(fix_i)-0.01, m_theta(fix_i),  m_theta(fix_i)+0.005, m_theta(fix_i)+0.01};     // Initial abscissae
	double *p_xinit = xinit;

	double xl = -2;
	double xr = 2;			  // Initial left and right (pseudo) extremes

    theta_params params;
    params.alpha = m_alpha;
    params.sum_failure = m_sum_failure_fix(fix_i);
    params.X_j = m_data->X.col(fix_i).cast<double>();  //Take from the fixed effects matrix
    params.epsilon = m_epsilon->array() + (params.X_j * m_theta(fix_i)).array(); // Adjust residual
    params.sigma_mu = m_sigma_mu;


	// Sample using ARS
    err = arms(xinit,ninit,&xl,&xr,theta_dens,&params,&convex,
			npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
	errorCheck(err);

    m_theta(fix_i) = xsamp[0];  // Save the new result
    *m_epsilon = params.epsilon - params.X_j * m_theta(fix_i); // Adjust residual
}

// Function for sampling marker effect (beta_i)
void BayesWBase::processColumn(const KernelPtr &kernel)
{
    auto * gaussKernel = dynamic_cast<BayesWKernel*>(kernel.get());
    assert(gaussKernel);

    const double beta_old = m_beta(gaussKernel->marker->i);

	//Change the residual vector only if the previous beta was non-zero
    if(beta_old != 0.0){
        *m_epsilon += *gaussKernel->calculateResidualUpdate(beta_old);
        //Also find the transformed residuals
        *m_vi = (m_alpha*m_epsilon->array()-EuMasc).exp();
	}

    gaussKernel->setVi(m_vi);
    gaussKernel->calculateSumFailure(m_failure_vector);

    /* Calculate the mixture probability */
    const double p = m_randomNumbers.at(kernel->marker->i);

    // Calculate the (ratios of) marginal likelihoods
    VectorXd marginal_likelihoods {m_K}; // likelihood for each mixture component
    // First element for the marginal likelihoods is always is pi_0 *sqrt(pi) for
    marginal_likelihoods(0) = m_pi_L(0) * sqrtPI;
    {
        const double exp_sum = gaussKernel->exponent_sum();

        for(int i=0; i < m_mixture_classes.size(); i++){
            //Calculate the sigma for the adaptive G-H
            double sigma = 1.0/sqrt(1 + m_alpha * m_alpha * m_sigma_b * m_mixture_classes(i) * exp_sum);
            marginal_likelihoods(i+1) = m_pi_L(i+1) * gauss_hermite_adaptive_integral(i, sigma, m_quad_points, gaussKernel);
        }
    }
	// Calculate the probability that marker is 0
    double acum = marginal_likelihoods(0)/marginal_likelihoods.sum();

    VectorXd localV = VectorXd::Zero(m_K);
    int component = 0;
    double beta_new = beta_old;
	//Loop through the possible mixture classes
    for (int k = 0; k < m_K; k++) {
		if (p <= acum) {
			//if zeroth component
			if (k == 0) {
                beta_new = 0;
			}
			// If is not 0th component then sample using ARS
			else {
                beta_params params;
                params.alpha = m_alpha;
                params.sigma_b = m_sigma_b;
                params.sum_failure = gaussKernel->sum_failure;
                params.used_mixture = m_mixture_classes(k-1);

                double safe_limit = 2 * sqrt(m_sigma_b * m_mixture_classes(k-1));

				// ARS parameters
                int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
                int neval;
                double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
                double convex = 1.0;
                int dometrop = 0;
                double xprev = 0.0;
                double xinit[4] = {beta_old - safe_limit/10 , beta_old,  beta_old + safe_limit/20, beta_old + safe_limit/10};     // Initial abscissae
                double *p_xinit = xinit;

                double xl = beta_old - safe_limit  ; //Construct the hull around previous beta value
                double xr = beta_old + safe_limit;

                // Sample using ARS
                err = estimateBeta(gaussKernel,m_epsilon,xinit,ninit,&xl,&xr, params, &convex,
                        npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
				errorCheck(err);

                beta_new = xsamp[0]; // Save the new result
			}

            localV[k] += 1.0;
            component = k;
			break;
		} else {
            if((k+1) == (m_K-1)){
				acum = 1; // In the end probability will be 1
			}else{
                acum += marginal_likelihoods(k+1)/marginal_likelihoods.sum();
			}
		}
	}

    // Only update m_epsilon if required
    const bool skipUpdate = beta_old == 0.0 && beta_new == 0.0;
    if (!skipUpdate) {
        //Re-update the residual vector
        *m_epsilon -= *gaussKernel->calculateResidualUpdate(beta_new);
        *m_vi = (m_alpha*m_epsilon->array()-EuMasc).exp();
    }

    m_v += localV;
    m_components(gaussKernel->marker->i) = component;
    m_beta(gaussKernel->marker->i) = beta_new;
}

// Function for sampling Weibull shape parameter (alpha)
void BayesWBase::sampleAlpha(){
	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
    double xinit[4] = {(m_alpha)*0.5, m_alpha,  (m_alpha)*1.15, (m_alpha)*1.5};     // Initial abscissae
	double *p_xinit = xinit;

	// Initial left and right (pseudo) extremes
	double xl = 0.0;
    double xr = 20.0;

    alpha_params params;
    params.alpha_0 = m_alpha_0;
    params.d = d;
    params.epsilon = *m_epsilon;
    params.failure_vector = m_failure_vector;
    params.kappa_0 = m_kappa_0;

	//Sample using ARS
    err = arms(xinit,ninit,&xl,&xr,alpha_dens,&params,&convex,
			npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
	errorCheck(err);
    m_alpha = xsamp[0];
}

/* Adaptive Gauss-Hermite version. Currently RAM solution */
int BayesWBase::runGibbs(AnalysisGraph* analysis)
{
    if (!analysis) {
        std::cout << "Cannot run BayesW analysis without a flow graph!" << std::endl;
        return 1;
    }

    const unsigned int M(m_data->numSnps);
    const unsigned int N(m_data->numInds);
    const unsigned int numFixedEffects(m_data->numFixedEffects);
    const int km1 = m_K - 1;

	init(M, N, numFixedEffects);
	int marker; //Marker index

	SampleWriter writer;
    writer.setFileName(m_outputFile);
	writer.setMarkerCount(M);
	writer.setIndividualCount(N);

	VectorXd sample(2*M+4); // variable containing a sample of all variables in the model: M marker effects, M mixture assignments, shape (alpha), mu, iteration number and sigma_b(sigma_g)

	//If we have fixed effects, then record their number to samplewriter and create a different header
	if(numFixedEffects > 0){
		writer.setFixedCount(numFixedEffects);
		writer.open_bayesW_fixed();
		sample.resize(numFixedEffects+2*M+4); // all the rest + theta (fixed effects)

	}else{
		writer.open_bayesW();
	}

	// Sampler variables
	std::vector<unsigned int> markerI(M);
	std::iota(markerI.begin(), markerI.end(), 0);

	std::cout<< "Running Gibbs sampling" << endl;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	// This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
	srand(2);
    for (int iteration = 0; iteration < m_max_iterations; iteration++) {

		/* 1. Intercept (mu) */
		sampleMu();

		/* 1a. Fixed effects (thetas) */
		if(numFixedEffects > 0){
			for(int fix_i = 0; fix_i < numFixedEffects; fix_i++){
				sampleTheta(fix_i);
			}
		}
		// Calculate the vector of exponent of the adjusted residuals
        *m_vi = (m_alpha*m_epsilon->array()-EuMasc).exp();

		std::random_shuffle(markerI.begin(), markerI.end());
		// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
		// 2. Sample beta parameters

		// Set counter for each mixture to be 1 ( (1,...,1) prior)
        m_v.setOnes();
        prepareForAnalysis();
        analysis->exec(this, N, M, markerI);

		// 3. Sample alpha parameter
		sampleAlpha();

		// 4. Sample sigma_b
        m_sigma_b = m_dist.inv_gamma_rng((double) (m_alpha_sigma + 0.5 * (M - m_v[0]+1)),
                (double)(m_beta_sigma + 0.5 * (M - m_v[0]+1) * m_beta.squaredNorm()));

		// 5. Sample prior mixture component probability from Dirichlet distribution
        m_pi_L = m_dist.dirichlet_rng(m_v.array());

		// Write the result to file
        if (iteration >= m_burn_in && iteration % m_thinning == 0) {
			if(numFixedEffects > 0){
                sample << iteration, m_alpha, m_mu, m_theta, m_beta,m_components.cast<double>(), m_sigma_b ;
			}else{
                sample << iteration, m_alpha, m_mu, m_beta,m_components.cast<double>(), m_sigma_b ;
			}
			writer.write(sample);
		}

		//Print results
        cout << iteration << ". " << M - m_v[0] +1 <<"; "<<m_v[1]-1 << "; "<<m_v[2]-1 << "; " << m_v[3]-1  <<"; " << m_alpha << "; " << m_sigma_b << endl;
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "duration: "<<duration << "s\n";

    return 0;
}

std::unique_ptr<AsyncResult> BayesWBase::processColumnAsync(const KernelPtr &kernel)
{
    auto * gaussKernel = dynamic_cast<BayesWKernel*>(kernel.get());
    assert(gaussKernel);

    // Local copies required to sample beta
    std::shared_ptr<VectorXd> epsilon;
    std::shared_ptr<VectorXd> vi;
    {
        std::shared_lock lock(m_mutex);
        epsilon = std::make_shared<VectorXd>(*m_epsilon);
        vi = std::make_shared<VectorXd>(*m_vi);
    }

    // No shared mutex for reading because no other thread writes to the values
    // specific to the marker this thread is working on

    const double beta_old = m_beta(gaussKernel->marker->i);

    auto result = std::make_unique<AsyncResult>();
    result->betaOld = beta_old;
    result->beta = beta_old;

    //Change the residual vector only if the previous beta was non-zero
    if(beta_old != 0.0){
        *epsilon += *gaussKernel->calculateResidualUpdate(beta_old);
        //Also find the transformed residuals
        *vi = (m_alpha*epsilon->array()-EuMasc).exp();
    }

    gaussKernel->setVi(vi);
    gaussKernel->calculateSumFailure(m_failure_vector);

    /* Calculate the mixture probability */
    const double p = m_randomNumbers.at(kernel->marker->i);

    // Calculate the (ratios of) marginal likelihoods
    VectorXd marginal_likelihoods {m_K}; // likelihood for each mixture component
    // First element for the marginal likelihoods is always is pi_0 *sqrt(pi) for
    marginal_likelihoods(0) = m_pi_L(0) * sqrtPI;
    {
        const double exp_sum = gaussKernel->exponent_sum();

        for(int i=0; i < m_mixture_classes.size(); i++){
            //Calculate the sigma for the adaptive G-H
            double sigma = 1.0/sqrt(1 + m_alpha * m_alpha * m_sigma_b * m_mixture_classes(i) * exp_sum);
            marginal_likelihoods(i+1) = m_pi_L(i+1) * gauss_hermite_adaptive_integral(i, sigma, m_quad_points, gaussKernel);
        }
    }
    // Calculate the probability that marker is 0
    double acum = marginal_likelihoods(0)/marginal_likelihoods.sum();

    result->v = std::make_unique<VectorXd>(VectorXd::Zero(m_K));
    int component = 0;
    //Loop through the possible mixture classes
    for (int k = 0; k < m_K; k++) {
        if (p <= acum) {
            //if zeroth component
            if (k == 0) {
                result->beta = 0;
            }
            // If is not 0th component then sample using ARS
            else {
                beta_params params;
                params.alpha = m_alpha;
                params.sigma_b = m_sigma_b;
                params.sum_failure = gaussKernel->sum_failure;
                params.used_mixture = m_mixture_classes(k-1);

                double safe_limit = 2 * sqrt(m_sigma_b * m_mixture_classes(k-1));

                // ARS parameters
                int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
                int neval;
                double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
                double convex = 1.0;
                int dometrop = 0;
                double xprev = 0.0;
                double xinit[4] = {beta_old - safe_limit/10 , beta_old,  beta_old + safe_limit/20, beta_old + safe_limit/10};     // Initial abscissae
                double *p_xinit = xinit;

                double xl = beta_old - safe_limit  ; //Construct the hull around previous beta value
                double xr = beta_old + safe_limit;

                // Sample using ARS
                err = estimateBeta(gaussKernel,epsilon,xinit,ninit,&xl,&xr, params, &convex,
                        npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
                errorCheck(err);

                result->beta = xsamp[0]; // Save the new result
            }

            (*result->v)(k) += 1.0;
            component = k;
            break;
        } else {
            if((k+1) == (m_K-1)){
                acum = 1; // In the end probability will be 1
            }else{
                acum += marginal_likelihoods(k+1)/marginal_likelihoods.sum();
            }
        }
    }

    // Only update m_epsilon if required
    const bool skipUpdate = result->betaOld == 0.0 && result->beta == 0.0;
    if (!skipUpdate) {
        result->deltaEpsilon = gaussKernel->calculateEpsilonChange(result->betaOld, result->beta);
    }

    m_components(gaussKernel->marker->i) = component;
    m_beta(gaussKernel->marker->i) = result->beta;

    return result;
}

void BayesWBase::doThreadSafeUpdates(const ConstAsyncResultPtr &result)
{
    assert(result);

    // No mutex required here - thread_safe_update_node is serial, therefore
    // only one runs at any time. m_v is not accessed elsewhere whilst the
    // flow graph is running.
    m_v += *result->v;
}

void BayesWBase::updateGlobal(const KernelPtr& kernel,
                              const ConstAsyncResultPtr& result)
{
    assert(kernel);
    assert(result);
    (void) kernel; // Unused

    std::unique_lock lock(m_mutex);

    *m_epsilon += *result->deltaEpsilon;
    *m_vi = (m_alpha*m_epsilon->array()-EuMasc).exp();
}
