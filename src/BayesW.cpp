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
#define sqrtPI 1.77245385090552
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

/* Function to check if ARS resulted with error*/
inline void errorCheck(int err){
	if(err>0){
		cout << "Error code = " << err << endl;
		exit(1);
	}
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

/* Function for the log density of some "fixed" covariate effect */
inline double theta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;
	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = - p.alpha * x * p.sum_failure - (((p.epsilon -  p.X_j * x)* p.alpha).array() - EuMasc).exp().sum() - x*x/(2*p.sigma_mu); // Prior is the same currently for intercepts and fixed effects
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

//The function for integration
inline double gh_integrand_adaptive(double s,double alpha, double dj, double sqrt_2Ck_sigmab, VectorXd vi, VectorXd Xj){
	//vi is a vector of exp(vi)
	double temp = -alpha *s*dj*sqrt_2Ck_sigmab + (vi.array()* (1 - (-Xj.array()*s*sqrt_2Ck_sigmab*alpha).exp() )).sum() -pow(s,2);
	return exp(temp);
}


//Calculate the value of the integral using Adaptive Gauss-Hermite quadrature
//Let's assume that mu is always 0 for speed
inline double gauss_hermite_adaptive_integral(int k, VectorXd vi,void *norm_data, double sigma, string n){
	pars p = *(static_cast<pars *>(norm_data));

	double temp = 0;
	double sqrt_2ck_sigma = sqrt(2*p.mixture_classes(k)*p.sigma_b);

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

		temp = 	w1 * gh_integrand_adaptive(x1,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w2 * gh_integrand_adaptive(x2,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
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

		temp = 	w1 * gh_integrand_adaptive(x1,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w2 * gh_integrand_adaptive(x2,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w3 * gh_integrand_adaptive(x3,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w4 * gh_integrand_adaptive(x4,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w5 ;//* gh_integrand_adaptive(x5,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j); // This part is just 1
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

		temp = 	w1 * gh_integrand_adaptive(x1,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w2 * gh_integrand_adaptive(x2,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w3 * gh_integrand_adaptive(x3,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w4 * gh_integrand_adaptive(x4,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w5 * gh_integrand_adaptive(x5,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w6 * gh_integrand_adaptive(x6,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w7;
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

		temp = 	w1 * gh_integrand_adaptive(x1,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w2 * gh_integrand_adaptive(x2,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w3 * gh_integrand_adaptive(x3,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w4 * gh_integrand_adaptive(x4,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w5 * gh_integrand_adaptive(x5,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w6 * gh_integrand_adaptive(x6,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w7 * gh_integrand_adaptive(x7,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w8 * gh_integrand_adaptive(x8,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w9 * gh_integrand_adaptive(x9,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w10 * gh_integrand_adaptive(x10,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j)+
				w11 ;//* gh_integrand_adaptive(x11,p.alpha,p.sum_failure,sqrt_2ck_sigma,vi,p.X_j);
	}else{
		cout << "Possible number of quad_points = 3,5,7,11" << endl;
		exit(1);
	}

	return sigma*temp;
}


//Pass the vector post_marginals of marginal likelihoods by reference
inline void marginal_likelihood_vec_calc(VectorXd prior_prob, VectorXd &post_marginals, VectorXd vi, void *norm_data, string n){
	pars p = *(static_cast<pars *>(norm_data));
	double exp_sum = (vi.array() * p.X_j.array() * p.X_j.array()).sum(); //For calculating sigma assume mu=0 and save time on computation
	// First element is pi_0 *sqrt(pi)
	post_marginals(0) = prior_prob(0) * sqrtPI;
	for(int i=0; i < p.mixture_classes.size(); i++){
		//Calculate the sigma for the adaptive G-H
		double sigma = 1.0/sqrt(1 + p.alpha * p.alpha * p.sigma_b * p.mixture_classes(i) * exp_sum);
		post_marginals(i+1) = prior_prob(i+1) * gauss_hermite_adaptive_integral(i, vi, norm_data, sigma, n);
	}
}


void BayesW::init(unsigned int markerCount, unsigned int individualCount, unsigned int fixedCount)
{
	// Read the failure indicator vector
	data.readFailureFile(opt.failureFile);

	// Component variables
	pi_L = VectorXd(K);           		 // prior mixture probabilities
	marginal_likelihoods = VectorXd(K);  // likelihood for each mixture component
	v = VectorXd(K);            		 // vector storing the component assignment

	// Linear model variables
	beta = VectorXd(markerCount);           // effect sizes
	theta = VectorXd(fixedCount);

	sum_failure = VectorXd(markerCount);	// Vector to sum SNP data vector * failure vector per SNP
	sum_failure_fix = VectorXd(fixedCount); // Vector to sum fixed vector * failure vector per fixed effect

	//phenotype vector
	y = VectorXd();

	// Resize the vectors in the structure
	used_data.X_j = VectorXd(individualCount);
	used_data.epsilon.resize(individualCount);
	used_data_alpha.epsilon.resize(individualCount);

	// Init the working variables
	const int km1 = K - 1;

	//vector with component class for each marker
	components = VectorXi(markerCount);
	components.setZero();

	//set priors for pi parameters
	//Give all mixtures (except 0 class) equal initial probabilities
	pi_L(0) = 0.99;
	pi_L.segment(1,km1).setConstant((1-pi_L(0))/km1);

	marginal_likelihoods.setOnes();   //Initialize with just ones

	beta.setZero();
	theta.setZero();

	//initialize epsilon vector as the phenotype vector
	y = data.y.cast<double>().array();
	mu = y.mean();       // mean or intercept

	// Initialize the variables in structures
	//Save variance classes
	used_data.mixture_classes.resize(km1);

	for(int i=0;i<(km1);i++){
		used_data.mixture_classes(i) = opt.S[i];   //Save the mixture data (C_k)
	}

	//Store the vector of failures only in the structure used for sampling alpha
	used_data_alpha.failure_vector = data.fail.cast<double>();

	double denominator = (6 * ((y.array() - mu).square()).sum()/(y.size()-1));
	used_data.alpha = PI/sqrt(denominator);    // The shape parameter initial value

	for(int i=0; i<(y.size()); ++i){
		(used_data.epsilon)[i] = y[i] - mu ; // Initially, all the BETA elements are set to 0, XBeta = 0
	}

	used_data.sigma_b = PI2/ (6 * pow(used_data.alpha,2) * markerCount ) ;

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

	// Save the number of events
	used_data.d = used_data_alpha.failure_vector.array().sum();
	used_data_alpha.d = used_data.d;

	// Save the sum(X_j*failure) for each j
	//Previous preprocessed version for reading columns
	//for(int marker=0; marker<M; marker++){
	//	sum_failure(marker) = ((data.mappedZ.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	//}
	for(int marker=0; marker < markerCount; marker++){
		sum_failure(marker) = ((data.Z.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
	}
	//If there are fixed effects, find the same values for them
	if(fixedCount > 0){
		for(int fix_i=0; fix_i < fixedCount; fix_i++){
			sum_failure_fix(fix_i) = ((data.X.col(fix_i).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
		}
	}
}
// Function for sampling intercept (mu)
void BayesW::sampleMu(){
	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
	double xinit[4] = {0.995*mu, mu,  1.005*mu, 1.01*mu};     // Initial abscissae
	double *p_xinit = xinit;

	double xl = 2;
	double xr = 5;   //xl and xr and the maximum and minimum values between which we sample

	used_data.epsilon = used_data.epsilon.array() + mu;// we add to epsilon =Y+mu-X*beta

	// Use ARS to sample mu (with density mu_dens, using parameters from used_data)
	err = arms(xinit,ninit,&xl,&xr,mu_dens,&used_data,&convex,
			npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

	errorCheck(err); // If there is error, stop the program
	mu = xsamp[0];   // Save the sampled value
	used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta
}

// Function for sampling fixed effect (theta_i)
void BayesW::sampleTheta(int fix_i){
	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
	double xinit[4] = {theta(fix_i)-0.01, theta(fix_i),  theta(fix_i)+0.005, theta(fix_i)+0.01};     // Initial abscissae
	double *p_xinit = xinit;

	double xl = -2;
	double xr = 2;			  // Initial left and right (pseudo) extremes

	used_data.X_j = data.X.col(fix_i).cast<double>();  //Take from the fixed effects matrix
	used_data.sum_failure = sum_failure_fix(fix_i);

	used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * theta(fix_i)).array(); // Adjust residual

	// Sample using ARS
	err = arms(xinit,ninit,&xl,&xr,theta_dens,&used_data,&convex,
			npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
	errorCheck(err);

	theta(fix_i) = xsamp[0];  // Save the new result
	used_data.epsilon = used_data.epsilon - used_data.X_j * theta(fix_i); // Adjust residual
}

// Function for sampling marker effect (beta_i)
void BayesW::sampleBeta(int marker){
	// Save the SNP effect column to a structure
	used_data.X_j = data.Z.col(marker).cast<double>();

	//Save sum(X_j*failure) to structure
	used_data.sum_failure = sum_failure(marker);

	//Change the residual vector only if the previous beta was non-zero
	if(beta(marker) != 0){
		used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
		//Also find the transformed residuals
		vi = (used_data.alpha*used_data.epsilon.array()-EuMasc).exp();
	}

	/* Calculate the mixture probability */
	double p = dist.unif_rng();  //Generate number from uniform distribution (for sampling from categorical distribution)

	// Calculate the (ratios of) marginal likelihoods
	marginal_likelihood_vec_calc(pi_L, marginal_likelihoods, vi, &used_data, quad_points);
	// Calculate the probability that marker is 0
	double acum = marginal_likelihoods(0)/marginal_likelihoods.sum();

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

				// ARS parameters
				int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
				int neval;
				double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
				double convex = 1.0;
				int dometrop = 0;
				double xprev = 0.0;
				double xinit[4] = {beta(marker) - safe_limit/10 , beta(marker),  beta(marker) + safe_limit/20, beta(marker) + safe_limit/10};     // Initial abscissae
				double *p_xinit = xinit;

				double xl = beta(marker) - safe_limit  ; //Construct the hull around previous beta value
				double xr = beta(marker) + safe_limit;

				// Sample using ARS
				err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
						npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
				errorCheck(err);

				beta(marker) = xsamp[0];  // Save the new result
				used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new
				vi = (used_data.alpha*used_data.epsilon.array()-EuMasc).exp();

				v[k] += 1.0;
				components[marker] = k;
			}
			break;
		} else {
			if((k+1) == (K-1)){
				acum = 1; // In the end probability will be 1
			}else{
				acum += marginal_likelihoods(k+1)/marginal_likelihoods.sum();
			}
		}
	}
}

// Function for sampling Weibull shape parameter (alpha)
void BayesW::sampleAlpha(){
	// ARS parameters
	int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
	int neval;
	double xsamp[0], xcent[10], qcent[10] = {5., 30., 70., 95.};
	double convex = 1.0;
	int dometrop = 0;
	double xprev = 0.0;
	double xinit[4] = {(used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3};     // Initial abscissae
	double *p_xinit = xinit;

	// Initial left and right (pseudo) extremes
	double xl = 0.0;
	double xr = 400.0;

	//Give the residual to alpha structure
	used_data_alpha.epsilon = used_data.epsilon;

	//Sample using ARS
	err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data_alpha,&convex,
			npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
	errorCheck(err);
	used_data.alpha = xsamp[0];
}

/* Adaptive Gauss-Hermite version. Currently RAM solution */
int BayesW::runGibbs_Gauss()
{
	const unsigned int M(data.numSnps);
	const unsigned int N(data.numInds);
	const unsigned int numFixedEffects(data.numFixedEffects);
	const int km1 = K - 1;

	init(M, N, numFixedEffects);
	int marker; //Marker index

	SampleWriter writer;
	writer.setFileName(outputFile);
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
	for (int iteration = 0; iteration < max_iterations; iteration++) {
		if (iteration > 0) {
			if (iteration % (int)std::ceil(max_iterations / 10) == 0)
				std::cout << "iteration: "<<iteration <<"\n";
		}
		/* 1. Intercept (mu) */
		sampleMu();

		/* 1a. Fixed effects (thetas) */
		if(numFixedEffects > 0){
			for(int fix_i = 0; fix_i < numFixedEffects; fix_i++){
				sampleTheta(fix_i);
			}
		}
		// Calculate the vector of exponent of the adjusted residuals
		vi = (used_data.alpha*used_data.epsilon.array()-EuMasc).exp();

		std::random_shuffle(markerI.begin(), markerI.end());
		// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
		// 2. Sample beta parameters

		// Set counter for each mixture to be 1 ( (1,...,1) prior)
		v.setOnes();
		for (int j = 0; j < M; j++) {
			marker = markerI[j];
			sampleBeta(marker);
		}
		// 3. Sample alpha parameter
		sampleAlpha();

		// 4. Sample sigma_b
		used_data.sigma_b = dist.inv_gamma_rng((double) (used_data.alpha_sigma + 0.5 * (M - v[0]+1)),
				(double)(used_data.beta_sigma + 0.5 * (M - v[0]+1) * beta.squaredNorm()));

		//Update the sqrt(2sigmab) variable
		used_data.sqrt_2sigmab = sqrt(2*used_data.sigma_b);

		// 5. Sample prior mixture component probability from Dirichlet distribution
		pi_L = dist.dirichilet_rng(v.array());

		// Write the result to file
		if ((iteration >= burn_in) and (iteration % thinning == 0)) {
			if(numFixedEffects > 0){
				sample << iteration, used_data.alpha, mu, theta, beta,components.cast<double>(), used_data.sigma_b ;
			}else{
				sample << iteration, used_data.alpha, mu, beta,components.cast<double>(), used_data.sigma_b ;
			}
			writer.write(sample);
		}

		//Print results
		cout << iteration << ". " << M - v[0] +1 <<"; "<<v[1]-1 << "; "<<v[2]-1 << "; " << v[3]-1  <<"; " << used_data.alpha << "; " << used_data.sigma_b << endl;
	}

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	std::cout << "duration: "<<duration << "s\n";

	return 0;
}


