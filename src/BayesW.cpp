/*
 * BayesW.cpp
 *
 *  Created on: 26 Nov 2018
 *  Author: Sven Erik Ojavee
 *  Last changes: 22 Feb 2019
 */

#include "data.hpp"
#include "distributions_boost.hpp"
//#include "concurrentqueue.h"
#include "options.hpp"
#include "BayesW.hpp"
#include "BayesW_arms.h"
#include "samplewriter.h"

#include <chrono>
#include <numeric>
#include <random>

/* Pre-calculate used constants */
#define PI 3.14159
#define PI2 6.283185
#define EuMasc 0.577215664901532

BayesW::BayesW(Data &data, Options &opt, const long memPageSize)
: seed(opt.seed)
, data(data)
, opt(opt)
, memPageSize(memPageSize)
, max_iterations(opt.chainLength)
, thinning(opt.thin)
, burn_in(opt.burnin)
, outputFile(opt.mcmcSampleFile)
, bedFile(opt.bedFile + ".bed")
, dist(opt.seed)
{

}


BayesW::~BayesW()
{
}



// Keep the necessary parameters in structures
// ARS uses the structure for using necessary parameters

struct pars{
	/* Common parameters for the densities */
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on
	//	 VectorXd epsilon_trunc;		// Difference of left truncation time and linear predictor (per subject)

	MatrixXd mixture_diff;    //Matrix to store (1/2Ck-1/2Cq) values
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

};

struct pars_alpha{
	VectorXd failure_vector;
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on

	/* Alpha-specific variables */
	double alpha_0, kappa_0;  /*  Prior parameters */

	/* Number of events (sum of failure indicators) */
	double d;
};


/* Function to assign initial values */
inline void assignArray(double *array_arg,VectorXd new_vals){
	for(int i = 0; i < new_vals.size(); i++){
		array_arg[i] = new_vals[i];
	}
}

/* Function to check if ARS resulted with error*/
inline void errorCheck(int err){
	if(err>0){
		cout << "Error code = " << err << endl;
		exit(1);
	}
}


/* Function for the second derivative. It is assumed that the residual (epsilon) is adjusted before */
/* C_k is the mixing proportion */
inline double beta_dens_der2(double x,double C_k, void *norm_data)
{
	double y;
	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = -(1/(C_k*p.sigma_b)) - pow(p.alpha,2) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp() ) *
			(p.X_j).array() * (p.X_j).array()).sum();
	return y;
};

/* Function for the ratio of der2 and der1 */
inline double beta_dens_12_ratio(double x, void *norm_data){

	pars p = *(static_cast<pars *>(norm_data));

	/*	return (-(x/(p.sigma_b)) - p.alpha * (p.X_j.array() * p.failure_vector.array()).sum() + (p.alpha)* ((((( p.epsilon * p.alpha).array() - (p.X_j * p.alpha).array() * x) - EuMasc).exp()) *
			(p.X_j).array()).sum())/
			(-(1/(p.sigma_b)) - (p.alpha)*(p.alpha) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp()) *
					(p.X_j).array() * (p.X_j).array()).sum());*/


	VectorXd exp_vector= (p.epsilon.array() - p.X_j .array() * x)* p.alpha - EuMasc;

	double max_val = exp_vector.maxCoeff();
	if(max_val > 700){
		// Subtract maximum from the vector (we are calculating the ratio assuming that , thus it does not change the final value)
		exp_vector = exp_vector.array() - max_val;
		exp_vector = exp_vector.array().exp();
		// Part with the failure vectors is assumed to be 0 now
		return((exp_vector.array() * p.X_j.array()).sum() / (exp_vector.array() * p.X_j.array() * p.X_j.array()).sum() / p.alpha);
	}

	exp_vector = exp_vector.array().exp();


	// This solution adds speed but the impact on power is not clear
	//	return (-(x/(p.sigma_b)) -  p.sum_failure + (exp_vector.array() * p.X_j.array()).sum())/
	//				(-(1/(p.sigma_b)) - (p.alpha)  * (exp_vector.array() * (p.X_j).array() * (p.X_j).array()).sum());

	return (-(  x/(p.sigma_b)) - p.alpha * p.sum_failure + p.alpha*(exp_vector.array() * p.X_j.array()).sum())/
			(-( 1/p.sigma_b) - (p.alpha)*(p.alpha)  * (exp_vector.array() * p.X_j.array() * p.X_j.array()).sum());

}

/* Function for the log density of mu */
inline double mu_dens(double x, void *norm_data)
/* We are sampling mu (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = - p.alpha * x * p.d - (( (p.epsilon).array()  - x) * p.alpha - EuMasc).exp().sum() - x*x/(2*p.sigma_mu);
	return y;
};

/* Function for the log density of alpha */
inline double alpha_dens(double x, void *norm_data)
/* We are sampling alpha (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars_alpha p = *(static_cast<pars_alpha *>(norm_data));
	y = (p.alpha_0 + p.d - 1) * log(x) + x * ((p.epsilon.array() * p.failure_vector.array()).sum() - p.kappa_0) -
			((p.epsilon * x).array() - EuMasc).exp().sum() ;
	return y;
};

/* Function for the log density of beta: uses mixture component from the structure norm_data */
inline double beta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	y = -p.alpha * x * p.sum_failure - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
			x * x / (2 * p.mixture_classes(p.used_mixture) * p.sigma_b) ;
	return y;
};

/* Function for the log density of beta: uses mixture component from C_k */
inline double beta_dens_ck(double x,double C_k, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	y = -p.alpha * x * p.sum_failure - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
			x * x / (2 * C_k * p.sigma_b) ;
	return y;
};

/* Function for the log density of beta evaluated at 0 */
inline double beta_dens_0(void *norm_data)
/* beta's log density evaluated at x=0*/
{
	double y;
	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y =  - ((p.epsilon * p.alpha).array() - EuMasc).exp().sum();
	return y;
};

/* Function to calculate the mode of the beta_j (assumed that C_k=inf) */
inline double betaMode(double initVal, void *my_data , double error = 0.000001, int max_count = 20){
	double x_i = initVal;
	double x_i1 = initVal + 0.00001;
	int counter = 0;

	while(abs(x_i-x_i1) > error){
		++counter;
		if(counter > max_count){
			return initVal;  //Failure if we repeat iterations too many times
		}
		x_i1 = x_i;
		x_i = x_i1 - beta_dens_12_ratio(x_i1,my_data);
	}
	return x_i;
}






///////////////////////////////////////////
/* Similar functions for left truncation */
///////////////////////////////////////////

/* Function for the log density of mu (LT)*/
/*
inline double mu_dens_ltrunc(double x, void *norm_data)
{
	double y;

	// In C++ we need to do a static cast for the void data
	pars p = *(static_cast<pars *>(norm_data));

	// cast voided pointer into pointer to struct norm_parm
	y = - p.alpha * x * p.failure_vector.sum() - (( (p.epsilon * p.alpha).array()  -  p.alpha * x) - EuMasc).exp().sum() +
			(( (p.epsilon_trunc * p.alpha).array()  -  p.alpha * x) - EuMasc).exp().sum() - x*x/(2*p.sigma_mu);
	return y;
};

// Function for the log density of alpha (LT)
inline double alpha_dens_ltrunc(double x, void *norm_data)
{
	double y;
	// In C++ we need to do a static cast for the void data
	pars p = *(static_cast<pars *>(norm_data));

	y = (p.alpha_0 + p.failure_vector.sum() - 1) * log(x) + x * ((p.epsilon.array() * p.failure_vector.array()).sum() - p.kappa_0) -
			((p.epsilon * x).array() - EuMasc).exp().sum() + ((p.epsilon_trunc * x).array() - EuMasc).exp().sum() ;
	return y;
};

// Function for the log density of beta (LT)
inline double beta_dens_ltrunc(double x, void *norm_data)
{
	double y;
	pars p = *(static_cast<pars *>(norm_data));

	y = -p.alpha * x * ((p.X_j).array() * (p.failure_vector).array()).sum() - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() +
			(((p.epsilon_trunc - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
			x * x / (2 * p.sigma_b) ;
	return y;
};

// Function calculates the difference beta_dens(0) - beta_dens(x)
inline double beta_dens_diff_ltrunc(double x, void *norm_data){
	double y;
	pars p = *(static_cast<pars *>(norm_data));
	y =   - ((p.epsilon * p.alpha).array() - EuMasc).exp().sum() +
			((p.epsilon_trunc * p.alpha).array() - EuMasc).exp().sum() -
			(-p.alpha * x * ((p.X_j).array() * (p.failure_vector).array()).sum() - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() +
					(((p.epsilon_trunc - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
					x * x / (2 * p.sigma_b) );

	return y;
}

// Function for the second derivative. It is assumed that the residual is adjusted before
inline double beta_dens_der2_ltrunc(double x, void *norm_data)
{
	double y;

	pars p = *(static_cast<pars *>(norm_data));

	y = -(1/(p.sigma_b)) + pow(p.alpha,2) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp() * (-1) +
			((( p.epsilon_trunc * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp()) *
			(p.X_j).array() * (p.X_j).array()).sum();
	return y;
};

// Function for the ratio of der2 and der1
inline double beta_dens_12_ratio_ltrunc(double x, void *norm_data){

	pars p = *(static_cast<pars *>(norm_data));

	return (-(x/(p.sigma_b)) - p.alpha * (p.X_j.array() * p.failure_vector.array()).sum() + (p.alpha)* ((((( p.epsilon * p.alpha).array() - (p.X_j * p.alpha).array() * x) - EuMasc).exp() -
			((( p.epsilon_trunc * p.alpha).array() - (p.X_j * p.alpha).array() * x) - EuMasc).exp()) *
			(p.X_j).array()).sum())/
			(-(1/(p.sigma_b)) + pow(p.alpha,2) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp() * (-1) +
					((( p.epsilon_trunc * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp()) *
					(p.X_j).array() * (p.X_j).array()).sum());
}

// Function for Beta mode
inline double betaMode_ltrunc(double initVal, void *my_data ,double error = 0.000001, int max_count = 20){
	double x_i = initVal;
	double x_i1 = initVal + 0.01;
	int counter = 0;

	while(abs(x_i-x_i1) > error){
		++counter;
		if(counter > max_count){
			return initVal;  //Failure
		}
		x_i1 = x_i;
		x_i = x_i1 - beta_dens_12_ratio_ltrunc(x_i1,my_data);
	}
	return x_i;
}
 */

//Integral for hk - h(0)
//Integrate from -0.5 to 0.5
inline double simson_integral_0(int quad_points_amount, int k, VectorXd vi, void *norm_data){
	pars p = *(static_cast<pars *>(norm_data));

	//VectorXd quad_points(quad_points_amount+1);
	double integral_value = 0;
	double quad_point = -0.5 ;
	integral_value = integral_value + exp(-p.alpha*quad_point*p.sum_failure + (vi.array().exp()*
					(1- (-p.alpha * p.X_j.array() * quad_point).exp())).sum() -
							quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b));

	for(int i=1; i < quad_points_amount; i++){
		quad_point = -0.5 + i*(0.5-(-0.5))/(quad_points_amount-1);
		if(i % 2 == 0 ){
			integral_value = integral_value + 2*exp(-p.alpha*quad_point*p.sum_failure + (vi.array().exp()*
							(1- (-p.alpha * p.X_j.array() * quad_point).exp())).sum() -
									quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b));
		}else{
			integral_value = integral_value + 4*exp(-p.alpha*quad_point*p.sum_failure + (vi.array().exp()*
										(1- (-p.alpha * p.X_j.array() * quad_point).exp())).sum() -
												quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b));
		}
	}
	quad_point = 0.5;
	integral_value = integral_value + exp(-p.alpha*quad_point*p.sum_failure + (vi.array().exp()*
						(1- (-p.alpha * p.X_j.array() * quad_point).exp())).sum() -
								quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b));

	return integral_value/(3*quad_points_amount);
}

inline double simson_integral_1(int quad_points_amount, int k, VectorXd vi, void *norm_data){
	pars p = *(static_cast<pars *>(norm_data));

	//VectorXd quad_points(quad_points_amount+1);
	double integral_value = 0;
	double quad_point = -0.5 ;
	integral_value = integral_value + exp(-p.alpha*quad_point*p.sum_failure -
			 (vi.array()-p.alpha * p.X_j.array() * quad_point).exp()).sum() -
					quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b);

	for(int i=1; i < quad_points_amount; i++){
		quad_point = -0.5 + i*(0.5-(-0.5))/(quad_points_amount-1);
		if(i % 2 == 0 ){
			integral_value = integral_value + 2*exp(-p.alpha*quad_point*p.sum_failure -
							 (vi.array()-p.alpha * p.X_j.array() * quad_point).exp()).sum() -
									quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b);
		}else{
			integral_value = integral_value + 4*exp(-p.alpha*quad_point*p.sum_failure -
					 (vi.array()-p.alpha * p.X_j.array() * quad_point).exp()).sum() -
							quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b);
		}
	}
	quad_point = 0.5;
	integral_value = integral_value + exp(-p.alpha*quad_point*p.sum_failure -
			 (vi.array()-p.alpha * p.X_j.array() * quad_point).exp()).sum() -
					quad_point*quad_point/(2*p.mixture_classes(k)*p.sigma_b);

	return integral_value/(3*quad_points_amount);
}

//Integrand part f(s) in the Gauss-Hermite quadrature formula
inline double gh_integrand(double s,double alpha, double dj, double sqrt_2Ck_sigmab, VectorXd vi, VectorXd Xj){
	double temp = -alpha *s*dj*sqrt_2Ck_sigmab - (vi.array() -  Xj.array()*s*sqrt_2Ck_sigmab*alpha ).exp().sum() +
			vi.array().exp().sum();
	return exp(temp);
}
//Same function evaluated at 0
/*inline double gh_integrand0(VectorXd vi){
	double temp = -(vi.array()).exp().sum();
	return exp(temp);
}*/

//Calculate numerically the value of marginal likelihood using Gauss-Hermite quadrature
//with n=5 points. Should be increased in the future
inline double gauss_hermite_integral(int k, VectorXd vi,void *norm_data){
	//At the moment specify all the quadrature points and weights manually
	//double x1 = 0;
	double w1 = 0.945309;

	double x2 = 0.958572;
	double w2 = 0.393619;
	double x3 = x2 * (-1);
	double w3 = w2;

	double x4 = 2.02018;
	double w4 = 0.0199532;
	double x5 = x4 * (-1);
	double w5 = w4;

	pars p = *(static_cast<pars *>(norm_data));

	double temp = w1 * 1 +  //gh_integrand(0) = 1
			w2 * gh_integrand(x2,p.alpha,p.sum_failure,sqrt(2*p.mixture_classes(k)*p.sigma_b),vi,p.X_j)+
			w3 * gh_integrand(x3,p.alpha,p.sum_failure,sqrt(2*p.mixture_classes(k)*p.sigma_b),vi,p.X_j)+
			w4 * gh_integrand(x4,p.alpha,p.sum_failure,sqrt(2*p.mixture_classes(k)*p.sigma_b),vi,p.X_j)+
			w5 * gh_integrand(x5,p.alpha,p.sum_failure,sqrt(2*p.mixture_classes(k)*p.sigma_b),vi,p.X_j);
	return temp;
}

inline double prob_calc0_gauss(VectorXd prior_prob, VectorXd vi, void *norm_data){
	double prob_0 = prior_prob(0) * sqrt(PI);

	pars p = *(static_cast<pars *>(norm_data));

	//Sum the comparisons
	for(int i=0; i < p.mixture_classes.size(); i++){
		prob_0 = prob_0 + prior_prob(i+1)* gauss_hermite_integral(i,vi,norm_data);
	}
	return prior_prob(0) * sqrt(PI)/prob_0;
}


inline double prob_calc0_simpson(int quad_points,VectorXd prior_prob, double C_0, VectorXd vi, void *norm_data){
	pars p = *(static_cast<pars *>(norm_data));
	double prob_0 = prior_prob(0) * C_0;

	//Sum the comparisons
	for(int i=0; i < p.mixture_classes.size(); i++){
		double simson_result = simson_integral_0(quad_points,i,vi,norm_data);
		prob_0 = prob_0 + prior_prob(i+1)*simson_result /sqrt(p.mixture_classes(i));
	}
	return prior_prob(0)*C_0/prob_0;
}

// Function that calculates probability of excluding the marker from the model */
inline double prob_calc0(double BETA_MODE, VectorXd prior_prob, double C_0, void *norm_data){

	pars p = *(static_cast<pars *>(norm_data));
	double prob_0 = prior_prob(0) * C_0;
	double beta_0 = beta_dens_0(norm_data);

	//Sum the comparisons
	for(int i=0; i < p.mixture_classes.size(); i++){
		prob_0 = prob_0 + prior_prob(i+1) * sqrt(-PI2/beta_dens_der2(BETA_MODE, p.mixture_classes(i), norm_data)/p.mixture_classes(i))*
				exp(beta_dens_ck(BETA_MODE,p.mixture_classes(i),norm_data)-beta_0);
	}
	return prior_prob(0)*C_0/prob_0;
}

/* Function that calculates probability of placing the marker into k-th mixture. Used if marker is included to the model */
inline double prob_calc(int k, double BETA_MODE, VectorXd prior_prob, double C_0, void *norm_data){
	pars p = *(static_cast<pars *>(norm_data));
	double beta_dens_der2_k = beta_dens_der2(BETA_MODE, p.mixture_classes(k), norm_data); //Calculate k-th second derivative
	double prob_k = prior_prob(0) * C_0 * sqrt(-beta_dens_der2_k/PI2)*exp(beta_dens_0(norm_data)-
			beta_dens_ck(BETA_MODE,p.mixture_classes(k),norm_data)) ;  //prior_prob vector has also the 0 component

	//Sum the comparisons
	for(int i=0; i<p.mixture_classes.size(); i++){
		prob_k = prob_k + prior_prob(i+1) * sqrt(beta_dens_der2_k/beta_dens_der2(BETA_MODE, p.mixture_classes(i), norm_data)/p.mixture_classes(i))*
				exp(pow(BETA_MODE,2)* p.mixture_diff(i,k)/p.sigma_b );  // We have previously calculated the differences to matrix
	}

	return prior_prob(k+1)/sqrt(p.mixture_classes(k))/prob_k;
}

inline double prob_calc0_marginal(VectorXd prior_prob, double U2, double V, void *norm_data){
	pars p = *(static_cast<pars *>(norm_data));
	double prob_0 = prior_prob(0)*sqrt(2*p.sigma_b);

	for(int i=0; i < p.mixture_classes.size(); i++){
			double V_i = V + 1/(2*p.sigma_b*p.mixture_classes(i));
			prob_0 = prob_0 + prior_prob(i+1)/sqrt(p.mixture_classes(i)*V_i)*
					exp(U2/(4*V_i)) ;
	}
	return prior_prob(0)*sqrt(2*p.sigma_b)/prob_0;
}

inline double prob_calc_marginal(int k, VectorXd prior_prob, double U2, double V, void *norm_data){
	pars p = *(static_cast<pars *>(norm_data));
	double prob_k = prior_prob(0)*sqrt(2*p.sigma_b);

	double k_likelihood;

	for(int i=0; i < p.mixture_classes.size(); i++){
		double V_i = V + 1/(2*p.sigma_b*p.mixture_classes(i));
		double added_amount = prior_prob(i+1)/sqrt(p.mixture_classes(i)*V_i)*exp(U2/(4*V_i));
		prob_k = prob_k + added_amount;
		if(i == k){  //Save the k-th element from the sum so we wouldn't have to calculate it again
			k_likelihood = added_amount;
		}
	}
	return k_likelihood/prob_k;
}

/* Functions to run each of the versions. Currently maintained one is runGibbs_Preprocessed */


/* Marginal Taylor version. Currently RAM solution */
int BayesW::runGibbs_Taylor()
{
	const unsigned int M(data.numSnps);
	const unsigned int N(data.numInds);
	const int K = opt.S.size()+1;  //number of mixtures + 0 class
	const int km1 = K -1;

	SampleWriter writer;
	writer.setFileName(outputFile);
	writer.setMarkerCount(M);
	writer.setIndividualCount(N);
	writer.open_bayesW();

	// Sampler variables
	VectorXd sample(2*M+5+K); // variable containing a sample of all variables in the model, M marker effects, shape (alpha), incl. prob (pi), mu, iteration number and beta variance,sigma_g
	std::vector<unsigned int> markerI(M);
	std::iota(markerI.begin(), markerI.end(), 0);

	data.readFailureFile(opt.failureFile);

	std::cout<< "Running Gibbs sampling" << endl;

	/* ARS parameters */
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};

	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;

	double xl, xr ;			  // Initial left and right (pseudo) extremes
	double xinit[4] = {2.5,3,5,10};     // Initial abscissae
	double *p_xinit = xinit;
	VectorXd new_xinit(4);  // Temporary vector to renew the initial parameters

	/* For ARS, we are keeping the data in this structure */
	struct pars used_data;
	struct pars_alpha used_data_alpha; // For alpha we keep it in a separate structure

	//mean and residual variables
	double mu;         // mean or intercept
	double BETA_MODE;  //Beta mode at hand


	//Save variance classes
	used_data.mixture_classes.resize(km1);

	for(int i=0;i<(km1);i++){
		used_data.mixture_classes(i) = opt.S[i];   //Save the mixture data (C_k)
	}
	// Component variables
	VectorXd pi_L(K); // Vector of mixture probabilities (+1 for 0 class)

	//Give all mixtures (except 0 class) equal initial probabilities
	pi_L(0) = 0.95;
	pi_L.segment(1,K-1).setConstant((1-pi_L(0))/km1);

	//Vector to contain probabilities of belonging to a mixture
	double acum;
	double betasqn = 0;

	VectorXd prob_vec(km1);   //exclude 0 mixture which is done before
	VectorXd v(K);            // variable storing the count in each component assignment

	//linear model variables   //y is logarithmed
	VectorXd beta(M);      // effect sizes

	int marker; //Marker index

	VectorXd gi(N); // The genetic effects vector
	gi.setZero();
	double sigma_g;

	VectorXd y;   //y is logarithmed here

	y = data.y.cast<double>();

	used_data_alpha.failure_vector = data.fail.cast<double>();

	beta.setZero(); //Exclude everything in the beginning

	//Initial value for intercept is the mean of the logarithms
	mu = y.mean();
	double denominator = (6 * ((y.array() - mu).square()).sum()/(y.size()-1));
	used_data.alpha = PI/sqrt(denominator);    // The shape parameter initial value

	(used_data.epsilon).resize(y.size());
	used_data_alpha.epsilon.resize(y.size());
	for(int i=0; i<(y.size()); ++i){
		(used_data.epsilon)[i] = y[i] - mu ; // Initially, all the BETA elements are set to 0, XBeta = 0
	}

	used_data.sigma_b = PI2/ (6 * pow(used_data.alpha,2) * M ) ;

	// Save the sum(X_j*failure) for each j
	VectorXd sum_failure(M);
	//for(int marker=0; marker<M; marker++){
	//	sum_failure(marker) = ((data.mappedZ.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	//}
	for(int marker=0; marker<M; marker++){
			sum_failure(marker) = ((data.Z.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	}

	// Save the number of events
	used_data.d = used_data_alpha.failure_vector.array().sum();
	used_data_alpha.d = used_data.d;

	/* Prior value selection for the variables */
	/* At the moment we set them to be weakly informative (in .hpp file) */
	/* alpha */
	used_data_alpha.alpha_0 = alpha_0;
	used_data_alpha.kappa_0 = kappa_0;
	/* mu */
	used_data.sigma_mu = sigma_mu;
	/* sigma_b */
	used_data.alpha_sigma = alpha_sigma;
	used_data.beta_sigma = beta_sigma;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	// This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
	srand(2);

	VectorXi components(M);
	components.setZero();  //Exclude all the markers from the model

	for (int iteration = 0; iteration < max_iterations; iteration++) {
		if (iteration > 0) {
			if (iteration % (int)std::ceil(max_iterations / 10) == 0)
				std::cout << "iteration: "<<iteration <<"\n";
		}
		/* 1. Mu */
		xl = 3; xr = 5;
		new_xinit << 0.95*mu, mu,  1.05*mu, 1.1*mu;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);
		used_data.epsilon = used_data.epsilon.array() + mu;//  we add the previous value

		err = arms(xinit,ninit,&xl,&xr,mu_dens,&used_data,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

		errorCheck(err);
		mu = xsamp[0];
		used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta

		std::random_shuffle(markerI.begin(), markerI.end());

		// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
		v.setOnes();           //Reset the counter
		double beta_diff_sum=0;
		for (int j = 0; j < M; j++) {
			marker = markerI[j];
			// Preprocessed solution
			//used_data.X_j = data.mappedZ.col(marker).cast<double>();

			// At the moment use RAM solution
			used_data.X_j = data.Z.col(marker).cast<double>();

			//Save sum(X_j*failure)
			used_data.sum_failure = sum_failure(marker);

			//Change the residual vector only if the previous beta was non-zero

			if(beta(marker) != 0){
				// Subtract the weighted last beta²
				used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
				betasqn = betasqn - beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
			}

			/* Calculate the mixture probability */
			double p = dist.unif_rng();  //Generate number from uniform distribution

			// Temporary variables for probability calculation
			VectorXd Xj_exp_vi = used_data.X_j.array() * (used_data.alpha*used_data.epsilon.array()-EuMasc).exp();
			double U2 = pow(used_data.alpha*(-used_data.sum_failure + Xj_exp_vi.sum()),2);
			double V = 0.5*pow(used_data.alpha,2) * (used_data.X_j.array() * Xj_exp_vi.array()).sum();

			// Calculate the probability that marker is 0
			acum = prob_calc0_marginal(pi_L, U2, V, &used_data);
			//Loop through the possible mixture classes
			for (int k = 0; k < K; k++) {
				if (p <= acum) {
					//if zeroth component
					if (k == 0) {
						beta(marker) = 0;
						v[k] += 1.0;
						components[marker] = k;
					}
					// If is not 0th component then sample using ARS
					else {
						used_data.used_mixture = k-1; // Save the mixture class before sampling (-1 because we count from 0)
						double safe_limit = 2 * sqrt(used_data.sigma_b * used_data.mixture_classes(k-1));
						xl = beta(marker) - safe_limit  ; //Construct the hull around previous beta value
						xr = beta(marker) + safe_limit;
						// Set initial values for constructing ARS hull
						new_xinit << beta(marker) - safe_limit/10 , beta(marker),  beta(marker) + safe_limit/20, beta(marker) + safe_limit/10;
						assignArray(p_xinit,new_xinit);
						// Sample using ARS
						err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
								npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
						errorCheck(err);

						beta(marker) = xsamp[0];  // Save the new result
						used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

						// Change the weighted sum of squares of betas
						v[k] += 1.0;
						components[marker] = k;
						betasqn = betasqn + beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
					}
					break;
				} else {
					if((k+1) == (K-1)){
						acum = 1; // In the end probability will be 1
					}else{
						acum += prob_calc_marginal(k, pi_L, U2, V, &used_data);

					}

				}
			}
		}
		// 3. Alpha
		xl = 0.0; xr = 400.0;
		new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);

		//Give the residual to alpha structure
		used_data_alpha.epsilon = used_data.epsilon;

		err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data_alpha,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
		errorCheck(err);
		used_data.alpha = xsamp[0];

		// 4. sigma_b
		if(true){
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
						(double)(used_data.beta_sigma + 0.5 * betasqn));
		}else{		//sigma_g version
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
			                            (double)(used_data.beta_sigma + 0.5 * (M - v[0]+1) * beta.squaredNorm()));
		}


		// 5. Mixture probability
		pi_L = dist.dirichilet_rng(v.array());
		// Also update the "spike parameter"

		if (iteration >= burn_in) {
			if (iteration % thinning == 0) {
				//6. Sigma_g
				//gi = y.array() - mu - used_data.epsilon.array();
				//sigma_g = (gi.array() * gi.array()).sum()/N - pow(gi.sum()/N,2);
				//sample << iteration, used_data.alpha, mu, beta,components, used_data.sigma_b , sigma_g;
				sample << iteration, used_data.alpha, mu, beta,components.cast<double>(), used_data.sigma_b ;
				writer.write(sample);
			}
		}

//Print results
        cout << iteration << ". " << M - v[0] +1 <<"; "<<v[1]-1 << "; "<<v[2]-1 << "; " << v[3]-1  <<"; " << used_data.alpha << "; " << used_data.sigma_b << endl;
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "duration: "<<duration << "s\n";

	return 0;
}



/* Marginal Taylor version. Currently RAM solution */
int BayesW::runGibbs_Gauss()
{
	const unsigned int M(data.numSnps);
	const unsigned int N(data.numInds);
	const int K = opt.S.size()+1;  //number of mixtures + 0 class
	const int km1 = K -1;

	SampleWriter writer;
	writer.setFileName(outputFile);
	writer.setMarkerCount(M);
	writer.setIndividualCount(N);
	writer.open_bayesW();

	// Sampler variables
	VectorXd sample(2*M+5+K); // variable containing a sample of all variables in the model, M marker effects, shape (alpha), incl. prob (pi), mu, iteration number and beta variance,sigma_g
	std::vector<unsigned int> markerI(M);
	std::iota(markerI.begin(), markerI.end(), 0);

	data.readFailureFile(opt.failureFile);

	std::cout<< "Running Gibbs sampling" << endl;

	/* ARS parameters */
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};

	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;

	double xl, xr ;			  // Initial left and right (pseudo) extremes
	double xinit[4] = {2.5,3,5,10};     // Initial abscissae
	double *p_xinit = xinit;
	VectorXd new_xinit(4);  // Temporary vector to renew the initial parameters

	/* For ARS, we are keeping the data in this structure */
	struct pars used_data;
	struct pars_alpha used_data_alpha; // For alpha we keep it in a separate structure

	//mean and residual variables
	double mu;         // mean or intercept
	double BETA_MODE;  //Beta mode at hand


	//Save variance classes
	used_data.mixture_classes.resize(km1);

	for(int i=0;i<(km1);i++){
		used_data.mixture_classes(i) = opt.S[i];   //Save the mixture data (C_k)
	}
	// Component variables
	VectorXd pi_L(K); // Vector of mixture probabilities (+1 for 0 class)

	//Give all mixtures (except 0 class) equal initial probabilities
	pi_L(0) = 0.95;
	pi_L.segment(1,K-1).setConstant((1-pi_L(0))/km1);

	//Vector to contain probabilities of belonging to a mixture
	double acum;
	double betasqn = 0;

	VectorXd prob_vec(km1);   //exclude 0 mixture which is done before
	VectorXd v(K);            // variable storing the count in each component assignment

	//linear model variables   //y is logarithmed
	VectorXd beta(M);      // effect sizes

	int marker; //Marker index

	VectorXd gi(N); // The genetic effects vector
	gi.setZero();
	double sigma_g;

	VectorXd y;   //y is logarithmed here

	y = data.y.cast<double>();

	used_data_alpha.failure_vector = data.fail.cast<double>();

	beta.setZero(); //Exclude everything in the beginning

	//Initial value for intercept is the mean of the logarithms
	mu = y.mean();
	double denominator = (6 * ((y.array() - mu).square()).sum()/(y.size()-1));
	used_data.alpha = PI/sqrt(denominator);    // The shape parameter initial value

	(used_data.epsilon).resize(y.size());
	used_data_alpha.epsilon.resize(y.size());
	for(int i=0; i<(y.size()); ++i){
		(used_data.epsilon)[i] = y[i] - mu ; // Initially, all the BETA elements are set to 0, XBeta = 0
	}

	used_data.sigma_b = PI2/ (6 * pow(used_data.alpha,2) * M ) ;

	// Save the sum(X_j*failure) for each j
	VectorXd sum_failure(M);
	//for(int marker=0; marker<M; marker++){
	//	sum_failure(marker) = ((data.mappedZ.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	//}
	for(int marker=0; marker<M; marker++){
			sum_failure(marker) = ((data.Z.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	}

	// Save the number of events
	used_data.d = used_data_alpha.failure_vector.array().sum();
	used_data_alpha.d = used_data.d;

	/* Prior value selection for the variables */
	/* At the moment we set them to be weakly informative (in .hpp file) */
	/* alpha */
	used_data_alpha.alpha_0 = alpha_0;
	used_data_alpha.kappa_0 = kappa_0;
	/* mu */
	used_data.sigma_mu = sigma_mu;
	/* sigma_b */
	used_data.alpha_sigma = alpha_sigma;
	used_data.beta_sigma = beta_sigma;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	// This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
	srand(2);

	VectorXi components(M);
	components.setZero();  //Exclude all the markers from the model

	for (int iteration = 0; iteration < max_iterations; iteration++) {
		if (iteration > 0) {
			if (iteration % (int)std::ceil(max_iterations / 10) == 0)
				std::cout << "iteration: "<<iteration <<"\n";
		}
		/* 1. Mu */
		xl = 3; xr = 5;
		new_xinit << 0.95*mu, mu,  1.05*mu, 1.1*mu;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);
		used_data.epsilon = used_data.epsilon.array() + mu;//  we add the previous value

		err = arms(xinit,ninit,&xl,&xr,mu_dens,&used_data,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

		errorCheck(err);
		mu = xsamp[0];
		used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta

		std::random_shuffle(markerI.begin(), markerI.end());

		// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
		v.setOnes();           //Reset the counter
		double beta_diff_sum=0;
		for (int j = 0; j < M; j++) {
			marker = markerI[j];
			// Preprocessed solution
			//used_data.X_j = data.mappedZ.col(marker).cast<double>();

			// At the moment use RAM solution
			used_data.X_j = data.Z.col(marker).cast<double>();

			//Save sum(X_j*failure)
			used_data.sum_failure = sum_failure(marker);

			//Change the residual vector only if the previous beta was non-zero

			if(beta(marker) != 0){
				// Subtract the weighted last beta²
				used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
				betasqn = betasqn - beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
			}

			/* Calculate the mixture probability */
			double p = dist.unif_rng();  //Generate number from uniform distribution

			// Temporary variables for probability calculation
			VectorXd vi = used_data.alpha*used_data.epsilon.array()-EuMasc;

			// Calculate the probability that marker is 0
			acum = prob_calc0_gauss(pi_L, vi, &used_data);




			//Loop through the possible mixture classes
			for (int k = 0; k < K; k++) {
				if (p <= acum) {
					//if zeroth component
					if (k == 0) {
						beta(marker) = 0;
						v[k] += 1.0;
						components[marker] = k;
					}
					// If is not 0th component then sample using ARS
					else {
						used_data.used_mixture = k-1; // Save the mixture class before sampling (-1 because we count from 0)
						double safe_limit = 2 * sqrt(used_data.sigma_b * used_data.mixture_classes(k-1));
						xl = beta(marker) - safe_limit  ; //Construct the hull around previous beta value
						xr = beta(marker) + safe_limit;
						// Set initial values for constructing ARS hull
						new_xinit << beta(marker) - safe_limit/10 , beta(marker),  beta(marker) + safe_limit/20, beta(marker) + safe_limit/10;
						assignArray(p_xinit,new_xinit);
						// Sample using ARS
						err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
								npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
						errorCheck(err);

						beta(marker) = xsamp[0];  // Save the new result
						used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

						// Change the weighted sum of squares of betas
						v[k] += 1.0;
						components[marker] = k;
						betasqn = betasqn + beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
					}
					break;
				} else {
					if((k+1) == (K-1)){
						acum = 1; // In the end probability will be 1
					}else{
						cout <<"Use only spike and slab --S 1" << endl;
						exit(1);

						//acum += prob_calc_marginal(k, pi_L, U2, V, &used_data);//Currently we try out only one mixture
					}

				}
			}
		}
		// 3. Alpha
		xl = 0.0; xr = 400.0;
		new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);

		//Give the residual to alpha structure
		used_data_alpha.epsilon = used_data.epsilon;

		err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data_alpha,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
		errorCheck(err);
		used_data.alpha = xsamp[0];

		// 4. sigma_b
		if(true){
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
						(double)(used_data.beta_sigma + 0.5 * betasqn));
		}else{		//sigma_g version
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
			                            (double)(used_data.beta_sigma + 0.5 * (M - v[0]+1) * beta.squaredNorm()));
		}


		// 5. Mixture probability
		pi_L = dist.dirichilet_rng(v.array());
		// Also update the "spike parameter"

		if (iteration >= burn_in) {
			if (iteration % thinning == 0) {
				//6. Sigma_g
				//gi = y.array() - mu - used_data.epsilon.array();
				//sigma_g = (gi.array() * gi.array()).sum()/N - pow(gi.sum()/N,2);
				//sample << iteration, used_data.alpha, mu, beta,components, used_data.sigma_b , sigma_g;
				sample << iteration, used_data.alpha, mu, beta,components.cast<double>(), used_data.sigma_b ;
				writer.write(sample);
			}
		}

//Print results
        cout << iteration << ". " << M - v[0] +1 <<"; "<<v[1]-1 << "; "<<v[2]-1 << "; " << v[3]-1  <<"; " << used_data.alpha << "; " << used_data.sigma_b << endl;
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "duration: "<<duration << "s\n";

	return 0;
}



/* Currently RAM solution */
int BayesW::runGibbs_old()
{
	const unsigned int M(data.numSnps);
	const unsigned int N(data.numInds);
	const int K = opt.S.size()+1;  //number of mixtures + 0 class
	const int km1 = K -1;

	//init();

	SampleWriter writer;
	writer.setFileName(outputFile);
	writer.setMarkerCount(M);
	writer.setIndividualCount(N);
	writer.open_bayesW();

	// Sampler variables
	VectorXd sample(2*M+5+K); // variable containing a sample of all variables in the model, M marker effects, shape (alpha), incl. prob (pi), mu, iteration number and beta variance,sigma_g
	std::vector<unsigned int> markerI(M);
	std::iota(markerI.begin(), markerI.end(), 0);

	data.readFailureFile(opt.failureFile);

	std::cout<< "Running Gibbs sampling" << endl;

	/* ARS parameters */
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};

	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;

	double xl, xr ;			  // Initial left and right (pseudo) extremes
	double xinit[4] = {2.5,3,5,10};     // Initial abscissae
	double *p_xinit = xinit;
	VectorXd new_xinit(4);  // Temporary vector to renew the initial parameters

	/* For ARS, we are keeping the data in this structure */
	struct pars used_data;
	struct pars_alpha used_data_alpha; // For alpha we keep it in a separate structure

	//mean and residual variables
	double mu;         // mean or intercept
	double BETA_MODE;  //Beta mode at hand

	//Precompute matrix of (1/2Ck - 1/2Cq)
	used_data.mixture_diff.resize(km1,km1);
	//Save variance classes
	used_data.mixture_classes.resize(km1);

	for(int i=0;i<(km1);i++){
		used_data.mixture_classes(i) = opt.S[i];   //Save the mixture data (C_k)
		for(int j=0;j<(km1);j++){
			used_data.mixture_diff(i,j) = 1/(2*opt.S[j]) - 1/(2*opt.S[i]);
		}
	}
	// Component variables
	VectorXd pi_L(K); // Vector of mixture probabilities (+1 for 0 class)
    VectorXd pi_L_cond1(km1); // Vector of conditional probabilities of belonging to specific mixture, given that the SNP has an effect

	//Give all mixtures (and 0 class) equal initial probabilities
	pi_L.setConstant(1.0/K);
    pi_L_cond1.setConstant(1.0/km1);

	//Vector to contain probabilities of belonging to a mixture
	double acum;
	double betasqn = 0;

	VectorXd prob_vec(km1);   //exclude 0 mixture which is done before
	VectorXd v(K);            // variable storing the count in each component assignment

	//linear model variables   //y is logarithmed
	VectorXd beta(M);      // effect sizes
	VectorXd BETAmodes(M); // Modes by variable
	BETAmodes.setZero();	//Initial value is 0

	int marker; //Marker index

	VectorXd gi(N); // The genetic effects vector
	gi.setZero();
	double sigma_g;

    double C_0 = (pi_L_cond1.array()/used_data.mixture_classes.array().sqrt()).array().sum() ; // Parameter for the estimating the belongingness to spike

	VectorXd y;   //y is logarithmed here

	y = data.y.cast<double>();

	used_data_alpha.failure_vector = data.fail.cast<double>();

	beta.setZero(); //Exclude everything in the beginning

	//Initial value for intercept is the mean of the logarithms
	mu = y.mean();
	double denominator = (6 * ((y.array() - mu).square()).sum()/(y.size()-1));
	used_data.alpha = PI/sqrt(denominator);    // The shape parameter initial value

	(used_data.epsilon).resize(y.size());
	used_data_alpha.epsilon.resize(y.size());
	for(int i=0; i<(y.size()); ++i){
		(used_data.epsilon)[i] = y[i] - mu ; // Initially, all the BETA elements are set to 0, XBeta = 0
	}

	used_data.sigma_b = PI2/ (6 * pow(used_data.alpha,2) * M ) ;

	// Save the sum(X_j*failure) for each j
	VectorXd sum_failure(M);
	//for(int marker=0; marker<M; marker++){
	//	sum_failure(marker) = ((data.mappedZ.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	//}
	for(int marker=0; marker<M; marker++){
			sum_failure(marker) = ((data.Z.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	}

	// Save the number of events
	used_data.d = used_data_alpha.failure_vector.array().sum();
	used_data_alpha.d = used_data.d;

	/* Prior value selection for the variables */
	/* At the moment we set them to be weakly informative (in .hpp file) */
	/* alpha */
	used_data_alpha.alpha_0 = alpha_0;
	used_data_alpha.kappa_0 = kappa_0;
	/* mu */
	used_data.sigma_mu = sigma_mu;
	/* sigma_b */
	used_data.alpha_sigma = alpha_sigma;
	used_data.beta_sigma = beta_sigma;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	// This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
	srand(2);

	VectorXi components(M);
	components.setZero();  //Exclude all the markers from the model

	for (int iteration = 0; iteration < max_iterations; iteration++) {
		if (iteration > 0) {
			if (iteration % (int)std::ceil(max_iterations / 10) == 0)
				std::cout << "iteration: "<<iteration <<"\n";
		}
		/* 1. Mu */
		xl = 3; xr = 5;
		new_xinit << 0.95*mu, mu,  1.05*mu, 1.1*mu;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);
		used_data.epsilon = used_data.epsilon.array() + mu;//  we add the previous value

		err = arms(xinit,ninit,&xl,&xr,mu_dens,&used_data,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

		errorCheck(err);
		mu = xsamp[0];
		used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta

		std::random_shuffle(markerI.begin(), markerI.end());

		// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
		v.setOnes();           //Reset the counter
		double beta_diff_sum=0;
		for (int j = 0; j < M; j++) {
			marker = markerI[j];
			// Preprocessed solution
			//used_data.X_j = data.mappedZ.col(marker).cast<double>();

			// At the moment use RAM solution
			used_data.X_j = data.Z.col(marker).cast<double>();

			//Save sum(X_j*failure)
			used_data.sum_failure = sum_failure(marker);

			//Change the residual vector only if the previous beta was non-zero

			if(beta(marker) != 0){
				// Subtract the weighted last beta²
				used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
				betasqn = betasqn - beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
			}

			/* Calculate the mixture probability */
			//VectorXd XjXj = used_data.X_j.array()*used_data.X_j.array();
			BETA_MODE = betaMode(BETAmodes(marker) ,&used_data);   //Find the posterior mode using the last mode as the starting value

			//TODO: For Multiple mixtures betaMode function needs to be defined
			BETAmodes(marker) = BETA_MODE;
			//TODO Calculate the second derivatives for each mixture component and save them

			double p = dist.unif_rng();  //Generate number from uniform distribution

			acum = prob_calc0(BETA_MODE,pi_L,C_0,&used_data);  // Calculate the probability that marker is 0
			//Loop through the possible mixture classes
			for (int k = 0; k < K; k++) {
				if (p <= acum) {
					//if zeroth component
					if (k == 0) {
						beta(marker) = 0;
						v[k] += 1.0;
						components[marker] = k;
					}
					// If is not 0th component then sample using ARS
					else {
						used_data.used_mixture = k-1; // Save the mixture class before sampling (-1 because we count from 0)
						double safe_limit = 2 * sqrt(used_data.sigma_b * used_data.mixture_classes(k-1));
						xl = BETA_MODE - safe_limit  ;
						xr = BETA_MODE + safe_limit;
						// Set initial values for constructing ARS hull
						new_xinit << BETA_MODE - safe_limit/10 , BETA_MODE,  BETA_MODE + safe_limit/20, BETA_MODE + safe_limit/10;
						assignArray(p_xinit,new_xinit);
						// Sample using ARS

						err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
								npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
						errorCheck(err);

						beta(marker) = xsamp[0];  // Save the result
						used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

						// Change the weighted sum of squares of betas
						v[k] += 1.0;
						components[marker] = k;
						betasqn = betasqn + beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
					}

					break;
				} else {
					acum += prob_calc(k,BETA_MODE,pi_L,C_0,&used_data);
				}
			}
		}
		// 3. Alpha
		xl = 0.0; xr = 400.0;
		new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);

		//Give the residual to alpha structure
		used_data_alpha.epsilon = used_data.epsilon;

		err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data_alpha,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
		errorCheck(err);
		used_data.alpha = xsamp[0];

		// 4. sigma_b
		if(true){
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
						(double)(used_data.beta_sigma + 0.5 * betasqn));
		}else{		//sigma_g version
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
			                            (double)(used_data.beta_sigma + 0.5 * (M - v[0]+1) * beta.squaredNorm()));
		}


		// 5. Mixture probability
		pi_L = dist.dirichilet_rng(v.array());
		// Also update the "spike parameter"
        pi_L_cond1 = pi_L.segment(1,km1).array()/(pi_L.segment(1,km1).array().sum());
        C_0 = (pi_L_cond1.array()/used_data.mixture_classes.array().sqrt()).array().sum() ;


		if (iteration >= burn_in) {
			if (iteration % thinning == 0) {
				//6. Sigma_g
				//gi = y.array() - mu - used_data.epsilon.array();
				//sigma_g = (gi.array() * gi.array()).sum()/N - pow(gi.sum()/N,2);
				//sample << iteration, used_data.alpha, mu, beta,components, used_data.sigma_b , sigma_g;
				sample << iteration, used_data.alpha, mu, beta,components.cast<double>(), used_data.sigma_b ;
				writer.write(sample);
			}
		}

//Print results
        cout << iteration << ". " << M - v[0] +1 <<"; "<<v[1]-1 << "; "<<v[2]-1 << "; " << v[3]-1  <<"; " << used_data.alpha << "; " << used_data.sigma_b << endl;
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "duration: "<<duration << "s\n";

	return 0;
}



/* Currently RAM solution */
//Use the idea from Laplace approximation but integrate with Simpson's method instead
int BayesW::runGibbs_Simpson()
{
	const unsigned int M(data.numSnps);
	const unsigned int N(data.numInds);
	const int K = opt.S.size()+1;  //number of mixtures + 0 class
	const int km1 = K -1;

	//init();

	SampleWriter writer;
	writer.setFileName(outputFile);
	writer.setMarkerCount(M);
	writer.setIndividualCount(N);
	writer.open_bayesW();

	// Sampler variables
	VectorXd sample(2*M+5+K); // variable containing a sample of all variables in the model, M marker effects, shape (alpha), incl. prob (pi), mu, iteration number and beta variance,sigma_g
	std::vector<unsigned int> markerI(M);
	std::iota(markerI.begin(), markerI.end(), 0);

	data.readFailureFile(opt.failureFile);

	std::cout<< "Running Gibbs sampling" << endl;

	/* ARS parameters */
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};

	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;

	double xl, xr ;			  // Initial left and right (pseudo) extremes
	double xinit[4] = {2.5,3,5,10};     // Initial abscissae
	double *p_xinit = xinit;
	VectorXd new_xinit(4);  // Temporary vector to renew the initial parameters

	/* For ARS, we are keeping the data in this structure */
	struct pars used_data;
	struct pars_alpha used_data_alpha; // For alpha we keep it in a separate structure

	//mean and residual variables
	double mu;         // mean or intercept

	//Precompute matrix of (1/2Ck - 1/2Cq)
	used_data.mixture_diff.resize(km1,km1);
	//Save variance classes
	used_data.mixture_classes.resize(km1);

	for(int i=0;i<(km1);i++){
		used_data.mixture_classes(i) = opt.S[i];   //Save the mixture data (C_k)
		for(int j=0;j<(km1);j++){
			used_data.mixture_diff(i,j) = 1/(2*opt.S[j]) - 1/(2*opt.S[i]);
		}
	}
	// Component variables
	VectorXd pi_L(K); // Vector of mixture probabilities (+1 for 0 class)
    VectorXd pi_L_cond1(km1); // Vector of conditional probabilities of belonging to specific mixture, given that the SNP has an effect

	//Give all mixtures (and 0 class) equal initial probabilities
	pi_L.setConstant(1.0/K);
    pi_L_cond1.setConstant(1.0/km1);

	//Vector to contain probabilities of belonging to a mixture
	double acum;
	double betasqn = 0;

	VectorXd prob_vec(km1);   //exclude 0 mixture which is done before
	VectorXd v(K);            // variable storing the count in each component assignment

	//linear model variables   //y is logarithmed
	VectorXd beta(M);      // effect sizes
	VectorXd BETAmodes(M); // Modes by variable
	BETAmodes.setZero();	//Initial value is 0

	int marker; //Marker index

	VectorXd gi(N); // The genetic effects vector
	gi.setZero();
	double sigma_g;

    double C_0 = (pi_L_cond1.array()/used_data.mixture_classes.array().sqrt()).array().sum() ; // Parameter for the estimating the belongingness to spike

	VectorXd y;   //y is logarithmed here

	y = data.y.cast<double>();

	used_data_alpha.failure_vector = data.fail.cast<double>();

	beta.setZero(); //Exclude everything in the beginning

	//Initial value for intercept is the mean of the logarithms
	mu = y.mean();
	double denominator = (6 * ((y.array() - mu).square()).sum()/(y.size()-1));
	used_data.alpha = PI/sqrt(denominator);    // The shape parameter initial value

	(used_data.epsilon).resize(y.size());
	used_data_alpha.epsilon.resize(y.size());
	for(int i=0; i<(y.size()); ++i){
		(used_data.epsilon)[i] = y[i] - mu ; // Initially, all the BETA elements are set to 0, XBeta = 0
	}

	used_data.sigma_b = PI2/ (6 * pow(used_data.alpha,2) * M ) ;

	// Save the sum(X_j*failure) for each j
	VectorXd sum_failure(M);

	for(int marker=0; marker<M; marker++){
			sum_failure(marker) = ((data.Z.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	}

	// Save the number of events
	used_data.d = used_data_alpha.failure_vector.array().sum();
	used_data_alpha.d = used_data.d;

	/* Prior value selection for the variables */
	/* At the moment we set them to be weakly informative (in .hpp file) */
	/* alpha */
	used_data_alpha.alpha_0 = alpha_0;
	used_data_alpha.kappa_0 = kappa_0;
	/* mu */
	used_data.sigma_mu = sigma_mu;
	/* sigma_b */
	used_data.alpha_sigma = alpha_sigma;
	used_data.beta_sigma = beta_sigma;

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	// This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
	srand(2);

	VectorXi components(M);
	components.setZero();  //Exclude all the markers from the model

	for (int iteration = 0; iteration < max_iterations; iteration++) {
		if (iteration > 0) {
			if (iteration % (int)std::ceil(max_iterations / 10) == 0)
				std::cout << "iteration: "<<iteration <<"\n";
		}
		/* 1. Mu */
		xl = 3; xr = 5;
		new_xinit << 0.95*mu, mu,  1.05*mu, 1.1*mu;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);
		used_data.epsilon = used_data.epsilon.array() + mu;//  we add the previous value

		err = arms(xinit,ninit,&xl,&xr,mu_dens,&used_data,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

		errorCheck(err);
		mu = xsamp[0];
		used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta

		std::random_shuffle(markerI.begin(), markerI.end());

		// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
		v.setOnes();           //Reset the counter
		double beta_diff_sum=0;
		for (int j = 0; j < M; j++) {
			marker = markerI[j];
			// Preprocessed solution
			//used_data.X_j = data.mappedZ.col(marker).cast<double>();

			// At the moment use RAM solution
			used_data.X_j = data.Z.col(marker).cast<double>();

			//Save sum(X_j*failure)
			used_data.sum_failure = sum_failure(marker);

			//Change the residual vector only if the previous beta was non-zero

			if(beta(marker) != 0){
				// Subtract the weighted last beta²
				used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
				betasqn = betasqn - beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
			}

			/* Calculate the mixture probability */
			double p = dist.unif_rng();  //Generate number from uniform distribution

			VectorXd vi = (used_data.alpha*used_data.epsilon.array()-EuMasc).exp();


			acum = prob_calc0_simpson(50,pi_L,C_0,vi,&used_data);  // Calculate the probability that marker is 0
			//Loop through the possible mixture classes
			for (int k = 0; k < K; k++) {
				if (p <= acum) {
					//if zeroth component
					if (k == 0) {
						beta(marker) = 0;
						v[k] += 1.0;
						components[marker] = k;
					}
					// If is not 0th component then sample using ARS
					else {
						used_data.used_mixture = k-1; // Save the mixture class before sampling (-1 because we count from 0)
						double safe_limit = 2 * sqrt(used_data.sigma_b * used_data.mixture_classes(k-1));
						xl = beta(marker) - safe_limit  ;
						xr = beta(marker) + safe_limit;
						// Set initial values for constructing ARS hull
						new_xinit << beta(marker) - safe_limit/10 , beta(marker),  beta(marker) + safe_limit/20, beta(marker) + safe_limit/10;
						assignArray(p_xinit,new_xinit);
						// Sample using ARS

						err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
								npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
						errorCheck(err);

						beta(marker) = xsamp[0];  // Save the result
						used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

						// Change the weighted sum of squares of betas
						v[k] += 1.0;
						components[marker] = k;
						betasqn = betasqn + beta(marker)*beta(marker)/used_data.mixture_classes(components[marker]-1);
					}

					break;
				} else {
					acum =1;//Right now only spike and slab
				}
			}
		}
		// 3. Alpha
		xl = 0.0; xr = 400.0;
		new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
		assignArray(p_xinit,new_xinit);

		//Give the residual to alpha structure
		used_data_alpha.epsilon = used_data.epsilon;

		err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data_alpha,&convex,
				npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
		errorCheck(err);
		used_data.alpha = xsamp[0];

		// 4. sigma_b
		if(true){
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
						(double)(used_data.beta_sigma + 0.5 * betasqn));
		}else{		//sigma_g version
			used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
			                            (double)(used_data.beta_sigma + 0.5 * (M - v[0]+1) * beta.squaredNorm()));
		}


		// 5. Mixture probability
		pi_L = dist.dirichilet_rng(v.array());
		// Also update the "spike parameter"
        pi_L_cond1 = pi_L.segment(1,km1).array()/(pi_L.segment(1,km1).array().sum());
        C_0 = (pi_L_cond1.array()/used_data.mixture_classes.array().sqrt()).array().sum() ;


		if (iteration >= burn_in) {
			if (iteration % thinning == 0) {
				//6. Sigma_g
				//gi = y.array() - mu - used_data.epsilon.array();
				//sigma_g = (gi.array() * gi.array()).sum()/N - pow(gi.sum()/N,2);
				//sample << iteration, used_data.alpha, mu, beta,components, used_data.sigma_b , sigma_g;
				sample << iteration, used_data.alpha, mu, beta,components.cast<double>(), used_data.sigma_b ;
				writer.write(sample);
			}
		}

//Print results
        cout << iteration << ". " << M - v[0] +1 <<"; "<<v[1]-1 << "; "<<v[2]-1 << "; " << v[3]-1  <<"; " << used_data.alpha << "; " << used_data.sigma_b << endl;
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "duration: "<<duration << "s\n";

	return 0;
}

/* PP solution with ability to handle left truncated data*/
/*
int BayesW::runGibbs_Preprocessed_LeftTruncated()
{
	int flag;
	moodycamel::ConcurrentQueue<Eigen::VectorXd> q;//lock-free queue
	const unsigned int M(data.numIncdSnps);
	const unsigned int N(data.numKeptInds);

	data.readFailureFile(opt.failureFile);
	data.readLeftTruncationFile(opt.leftTruncFile);


	VectorXi gamma(M);
	VectorXf normedSnpData(data.numKeptInds);

	flag = 0;

	std::cout<<"Running Gibbs sampling" << endl;

	// Compute the SNP data length in bytes
	size_t snpLenByt = (data.numInds % 4) ? data.numInds / 4 + 1 : data.numInds / 4;

	omp_set_nested(1); // 1 - enables nested parallelism; 0 - disables nested parallelism.

	//Eigen::initParallel();

#pragma omp parallel shared(flag,q)
	{
#pragma omp sections
		{

			{
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
				VectorXd new_xinit(4);  // Temporary vector to renew the initial parameters

				struct pars used_data;  // For ARS, we are keeping the data in this structure

				//mean and residual variables
				double mu; // mean or intercept

				double prob;  //Inclusion probability
				double BETA_MODE;  //Beta mode at hand

				//component variables
				double pi = 0.5; // prior inclusion probability

				//linear model variables   //y is logarithmed
				VectorXd beta(M); // effect sizes
				VectorXd BETAmodes(M); // Modes by variable

				//     VectorXd y_tilde(N); // variable containing the adjusted residual to exclude the effects of a given marker

				//sampler variables
				VectorXd sample(1 * M + 5); // variable containing a sample of all variables in the model, M marker effects, shape (alpha), incl. prob (pi), mu, iteration number and beta variance
				std::vector<int> markerI(M);
				std::iota(markerI.begin(), markerI.end(), 0);

				int marker;

				VectorXd y;   //y is logarithmed

				y = data.y.cast<double>();

				VectorXi failure(N);   //Failure vector
				failure = data.fail;

				VectorXd left_trunc(N);   //Left truncation time vector
				left_trunc = data.left_trunc.cast<double>();


				(used_data.epsilon).resize(y.size());
				(used_data.failure_vector).resize(failure.size());

				//         y_tilde.setZero();
				beta.setZero();
				//Initial value for intercept is the mean of the logarithms

				mu = y.mean();
				double denominator = (6 * ((y.array() - mu).square()).sum()/(y.size()-1));
				used_data.alpha = PI/sqrt(denominator);    // The shape parameter initial value


				gamma.setZero();  //Exclude all the markers from the model

				for(int i=0; i<(y.size()); ++i){
					(used_data.epsilon)[i] = y[i] - mu ; // Initially, all the BETA elements are set to 0, XBeta = 0
					(used_data.epsilon_trunc)[i] = left_trunc[i] - mu;  //Also use the "left truncation residual"
					(used_data.failure_vector)[i] = failure[i];
				}


				used_data.sigma_b = pow(PI,2)/ (6 * pow(used_data.alpha,2) * M ) ;


				// alpha
				used_data.alpha_0 = alpha_0;
				used_data.kappa_0 = kappa_0;
				// mu
				used_data.sigma_mu = sigma_mu;
				// sigma_b
				used_data.alpha_sigma = alpha_sigma;
				used_data.beta_sigma = beta_sigma;


				std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

				// Need to think whether log survival data should be scaled

				//             y = (data.y.cast<double>().array() - data.y.cast<double>().mean());
				//             y /= sqrt(y.squaredNorm() / ((double)N - 1.0));


				// This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
				srand(2);

				for (int iteration = 0; iteration < max_iterations; iteration++) {

					if (iteration > 0) {
						if (iteration % (int)std::ceil(max_iterations / 10) == 0)
							std::cout << "iteration: "<<iteration <<"\n";
					}



					// 1. Mu
					xl = -5; xr = 10.0;
					new_xinit << 0.95*mu, mu,  1.05*mu, 1.1*mu;  // New values for abscissae evaluation
					assignArray(p_xinit,new_xinit);
					used_data.epsilon = used_data.epsilon.array() + mu;//  we add the previous value
					used_data.epsilon_trunc = used_data.epsilon_trunc.array() + mu;

					err = arms(xinit,ninit,&xl,&xr,mu_dens_ltrunc,&used_data,&convex,
							npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

					errorCheck(err);
					mu = xsamp[0];
					used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta
					used_data.epsilon_trunc = used_data.epsilon_trunc.array() - mu;


					std::random_shuffle(markerI.begin(), markerI.end());

					// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
					for (int j = 0; j < M; j++) {

						marker = markerI[j];

						// Using the preprocessed solution
						used_data.X_j = data.mappedZ.col(marker).cast<double>();

						used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
						used_data.epsilon_trunc = used_data.epsilon_trunc.array() + (used_data.X_j * beta(marker)).array();


						// Calculate the inclusion probability
						if( true or (iteration <= burn_in)){ // currently use Newton's method always
							BETA_MODE = betaMode(beta(marker),&used_data);
							BETAmodes(marker) = BETA_MODE;
						}else{
							BETA_MODE = BETAmodes(marker);
						}

						prob = 1/(1 + ((1-pi)/pi) * exp(beta_dens_diff_ltrunc(BETA_MODE,&used_data)) *
								sqrt(-beta_dens_der2_ltrunc(BETA_MODE,&used_data)/(2*PI)));

						gamma(marker) = dist.bernoulli(prob);   // Sample the inclusion based on the probability

						if(gamma(marker) == 1){
							new_xinit << BETA_MODE - 0.1 , BETA_MODE,  BETA_MODE+0.05, BETA_MODE + 0.1;
							assignArray(p_xinit,new_xinit);
							err = arms(xinit,ninit,&xl,&xr,beta_dens_ltrunc,&used_data,&convex,
									npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
							errorCheck(err);
							beta(marker) = xsamp[0];
							used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new
							used_data.epsilon_trunc = used_data.epsilon_trunc - used_data.X_j * beta(marker);

						}else{
							beta(marker) = 0;
							//If beta is 0, then we don't need to do the residual updates anymore
						}


					}

					// 3. Alpha
					xl = 0.0; xr = 400.0;
					new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
					assignArray(p_xinit,new_xinit);


					err = arms(xinit,ninit,&xl,&xr,alpha_dens_ltrunc,&used_data,&convex,
							npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
					errorCheck(err);
					used_data.alpha = xsamp[0];

					// 4. sigma_b
					used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (gamma.sum())),(double)(used_data.beta_sigma + 0.5 * (beta.array() * beta.array()).sum()));

					//5. Inclusion probability
					pi = dist.beta_rng(1+gamma.sum(), 1 + gamma.size() - gamma.sum());

					if (iteration >= burn_in) {
						if (iteration % thinning == 0) {
							sample << iteration, used_data.alpha, mu, beta, used_data.sigma_b ,pi;
							q.enqueue(sample);
						}
					}

					cout << iteration << ". " << gamma.sum() <<"; " << used_data.alpha << endl;
				}

				std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
				std::cout << "duration: "<<duration << "s\n";
				flag = 1;

			}

			//this thread saves in the output file using the lock-free queue
#pragma omp section
			{
				bool queueFull;
				queueFull = 0;
				std::ofstream outFile;
				outFile.open(outputFile);
				VectorXd sampleq(1 * M + 5 );
				IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "");
				outFile<< "iteration," << "alpha," << "mu,";
				for (unsigned int i = 0; i < M; ++i) {
					outFile << "beta[" << (i+1) << "],";
				}
				outFile << "sigma_b," << "pi,";

				outFile << "\n";

				while (!flag) {
					if (q.try_dequeue(sampleq))
						outFile << sampleq.transpose().format(CommaInitFmt) << "\n";
				}
			}
		}
	}

	return 0;
}
 */
