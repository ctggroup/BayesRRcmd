/*
 * BayesW.cpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#include "data.hpp"
#include "distributions_boost.hpp"
#include "concurrentqueue.h"
#include "options.hpp"
#include "BayesW.hpp"
#include "BayesW_arms.h"


#include <numeric>
#include <random>

#define PI 3.14159
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
	float* ptr =(float*)&opt.S[0];
	//   cva=(Eigen::Map<Eigen::VectorXf>(ptr,opt.S.size())).cast<double>();
}


BayesW::~BayesW()
{
}

// Keep the necessary parameters in a structure
struct pars{
	/* Common parameters for the densities */
	VectorXd failure_vector;	// Data for failures (per subject)
	VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on
	VectorXd epsilon_trunc;		// Difference of left truncation time and linear predictor (per subject)

	/* Store the current variables */
	double alpha, sigma_b;

	/* Alpha-specific variables */
	double alpha_0, kappa_0;  /*  Prior parameters */
	/* Beta_j - specific variables */
	VectorXd X_j;
	/* Mu-specific variables */
	double sigma_mu;
	/* sigma_b-specific variables */
	double alpha_sigma, beta_sigma;

};

/* Function to assign initial values */
void assignArray(double *array_arg,VectorXd new_vals){
	for(int i = 0; i < new_vals.size(); i++){
		array_arg[i] = new_vals[i];
	}
}

void errorCheck(int err){
	if(err>0){
		cout << "Error code = " << err << endl;
		exit(1);
	}
}



double beta_dens_der2(double x, void *norm_data)
/* Function for the second derivative. It is assumed that the residual is adjusted before */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = -(1/(p.sigma_b)) - pow(p.alpha,2) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp() ) *
			(p.X_j).array() * (p.X_j).array()).sum();
	return y;
};

/* Function for the ratio of der2 and der1 */
double beta_dens_12_ratio(double x, void *norm_data){

	pars p = *(static_cast<pars *>(norm_data));

	return (-(x/(p.sigma_b)) - p.alpha * (p.X_j.array() * p.failure_vector.array()).sum() + (p.alpha)* ((((( p.epsilon * p.alpha).array() - (p.X_j * p.alpha).array() * x) - EuMasc).exp()) *
			(p.X_j).array()).sum())/
			(-(1/(p.sigma_b)) - (p.alpha)*(p.alpha) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp()) *
					(p.X_j).array() * (p.X_j).array()).sum());
}

/* Function calculates the difference beta_dens(0) - beta_dens(x) */
double beta_dens_diff(double x, void *norm_data){
	double y;
	pars p = *(static_cast<pars *>(norm_data));
	y =  - (((p.epsilon) * p.alpha).array() - EuMasc).exp().sum() -
			(-p.alpha * x * ((p.X_j).array() * (p.failure_vector).array()).sum() - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
					x * x / (2 * p.sigma_b) );

	return y;
}

double mu_dens(double x, void *norm_data)
/* mu's log density */
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = - p.alpha * x * p.failure_vector.sum() - (( (p.epsilon * p.alpha).array()  -  p.alpha * x) - EuMasc).exp().sum() - x*x/(2*p.sigma_mu);
	return y;
};

double alpha_dens(double x, void *norm_data)
/* alpha's log density */
/* We are sampling alpha (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */



	y = (p.alpha_0 + p.failure_vector.sum() - 1) * log(x) + x * ((p.epsilon.array() * p.failure_vector.array()).sum() - p.kappa_0) -
			((p.epsilon * x).array() - EuMasc).exp().sum() ;
	return y;
};

double beta_dens(double x, void *norm_data)
/* beta's log density */
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = -p.alpha * x * ((p.X_j).array() * (p.failure_vector).array()).sum() - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
			x * x / (2 * p.sigma_b) ;
	return y;
};

double betaMode(double initVal, void *my_data ,double error = 0.000001, int max_count = 20){
	double x_i = initVal;
	double x_i1 = initVal + 0.01;
	int counter = 0;

	while(abs(x_i-x_i1) > error){
		++counter;
		if(counter > max_count){
			return initVal;  //Failure
		}
		x_i1 = x_i;
		//x_i = x_i1 - beta_dens_der1(x_i1,my_data)/beta_dens_der2(x_i1,my_data);
		x_i = x_i1 - beta_dens_12_ratio(x_i1,my_data);
	}
	return x_i;
}


/* Similar functions for left truncation */
double mu_dens_ltrunc(double x, void *norm_data)
/* mu's log density */
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = - p.alpha * x * p.failure_vector.sum() - (( (p.epsilon * p.alpha).array()  -  p.alpha * x) - EuMasc).exp().sum() +
			 (( (p.epsilon_trunc * p.alpha).array()  -  p.alpha * x) - EuMasc).exp().sum() - x*x/(2*p.sigma_mu);
	return y;
};

double alpha_dens_ltrunc(double x, void *norm_data)
/* alpha's log density */
/* We are sampling alpha (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */


	y = (p.alpha_0 + p.failure_vector.sum() - 1) * log(x) + x * ((p.epsilon.array() * p.failure_vector.array()).sum() - p.kappa_0) -
			((p.epsilon * x).array() - EuMasc).exp().sum() + ((p.epsilon_trunc * x).array() - EuMasc).exp().sum() ;
	return y;
};

double beta_dens_ltrunc(double x, void *norm_data)
/* beta's log density */
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = -p.alpha * x * ((p.X_j).array() * (p.failure_vector).array()).sum() - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() +
			(((p.epsilon_trunc - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
			x * x / (2 * p.sigma_b) ;
	return y;
};
/* Function calculates the difference beta_dens(0) - beta_dens(x) */
double beta_dens_diff_ltrunc(double x, void *norm_data){
	double y;
	pars p = *(static_cast<pars *>(norm_data));
	y =   - ((p.epsilon * p.alpha).array() - EuMasc).exp().sum() +
			((p.epsilon_trunc * p.alpha).array() - EuMasc).exp().sum() -
			(-p.alpha * x * ((p.X_j).array() * (p.failure_vector).array()).sum() - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() +
					(((p.epsilon_trunc - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
					x * x / (2 * p.sigma_b) );

	return y;
}

double beta_dens_der2_ltrunc(double x, void *norm_data)
/* Function for the second derivative. It is assumed that the residual is adjusted before */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
	pars p = *(static_cast<pars *>(norm_data));

	/* cast voided pointer into pointer to struct norm_parm */
	y = -(1/(p.sigma_b)) + pow(p.alpha,2) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp() * (-1) +
			((( p.epsilon_trunc * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp()) *
			(p.X_j).array() * (p.X_j).array()).sum();
	return y;
};

/* Function for the ratio of der2 and der1 */
double beta_dens_12_ratio_ltrunc(double x, void *norm_data){

	pars p = *(static_cast<pars *>(norm_data));

	return (-(x/(p.sigma_b)) - p.alpha * (p.X_j.array() * p.failure_vector.array()).sum() + (p.alpha)* ((((( p.epsilon * p.alpha).array() - (p.X_j * p.alpha).array() * x) - EuMasc).exp() -
			((( p.epsilon_trunc * p.alpha).array() - (p.X_j * p.alpha).array() * x) - EuMasc).exp()) *
			(p.X_j).array()).sum())/
			(-(1/(p.sigma_b)) + pow(p.alpha,2) *  ((((( p.epsilon * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp() * (-1) +
					((( p.epsilon_trunc * p.alpha).array() - (p.X_j* p.alpha).array() * x) - EuMasc).exp()) *
					(p.X_j).array() * (p.X_j).array()).sum());
}

/* Function for Beta mode */
double betaMode_ltrunc(double initVal, void *my_data ,double error = 0.000001, int max_count = 20){
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


int BayesW::runGibbs_notPreprocessed()
{
	int flag;
	moodycamel::ConcurrentQueue<Eigen::VectorXd> q;//lock-free queue
	const unsigned int M(data.numIncdSnps);
	const unsigned int N(data.numKeptInds);

	data.readFailureFile(opt.failureFile);


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
					(used_data.failure_vector)[i] = failure[i];
				}


				used_data.sigma_b = pow(PI,2)/ (6 * pow(used_data.alpha,2) * M ) ;

				/* Prior value selection for the variables */
				/* At the moment we set them to be uninformative */
				/* alpha */
				used_data.alpha_0 = alpha_0;
				used_data.kappa_0 = kappa_0;
				/* mu */
				used_data.sigma_mu = sigma_mu;
				/* sigma_b */
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



					/* 1. Mu */
					xl = -5; xr = 10.0;
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
					for (int j = 0; j < M; j++) {

						marker = markerI[j];

						// For not preprocessed we do the following
						data.getSnpDataFromBedFileUsingMmap_openmp(bedFile, snpLenByt, memPageSize, marker, normedSnpData);
						//I use a temporal variable to do the cast, there should be better ways to do this.
						used_data.X_j = normedSnpData.cast<double>();


						used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();


						/* Calculate the inclusion probability */
						if( true or (iteration <= burn_in)){ //or some other variable should be used
							BETA_MODE = betaMode(beta(marker),&used_data);
							BETAmodes(marker) = BETA_MODE;
						}else{
							BETA_MODE = BETAmodes(marker);
						}

						prob = 1/(1 + ((1-pi)/pi) * exp(beta_dens_diff(BETA_MODE,&used_data)) *
								sqrt(-beta_dens_der2(BETA_MODE,&used_data)/(2*PI)));
						//cout << marker << ": "<< prob << ";";

						gamma(marker) = dist.bernoulli(prob);   // Sample the inclusion based on the probability

						if(gamma(marker) == 1){
							new_xinit << BETA_MODE - 0.1 , BETA_MODE,  BETA_MODE+0.05, BETA_MODE + 0.1;
							assignArray(p_xinit,new_xinit);
							err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
									npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
							errorCheck(err);
							beta(marker) = xsamp[0];
							used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

						}else{
							beta(marker) = 0;
							//If beta is 0, then we don't need to do the residual update anymore
						}


					}

					// 3. Alpha
					xl = 0.0; xr = 400.0;
					new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
					assignArray(p_xinit,new_xinit);


					err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data,&convex,
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


int BayesW::runGibbs_Preprocessed()
{
	int flag;
	moodycamel::ConcurrentQueue<Eigen::VectorXd> q;//lock-free queue
	const unsigned int M(data.numIncdSnps);
	const unsigned int N(data.numKeptInds);

	data.readFailureFile(opt.failureFile);


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
					(used_data.failure_vector)[i] = failure[i];
				}


				used_data.sigma_b = pow(PI,2)/ (6 * pow(used_data.alpha,2) * M ) ;

				/* Prior value selection for the variables */
				/* At the moment we set them to be uninformative */
				/* alpha */
				used_data.alpha_0 = alpha_0;
				used_data.kappa_0 = kappa_0;
				/* mu */
				used_data.sigma_mu = sigma_mu;
				/* sigma_b */
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



					/* 1. Mu */
					xl = -5; xr = 10.0;
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
					for (int j = 0; j < M; j++) {

						marker = markerI[j];

						// Using the preprocessed solution
						used_data.X_j = data.mappedZ.col(marker).cast<double>();

						used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();


						/* Calculate the inclusion probability */
						if( true or (iteration <= burn_in)){ //or some other variable should be used
							BETA_MODE = betaMode(beta(marker),&used_data);
							BETAmodes(marker) = BETA_MODE;
						}else{
							BETA_MODE = BETAmodes(marker);
						}

						prob = 1/(1 + ((1-pi)/pi) * exp(beta_dens_diff(BETA_MODE,&used_data)) *
								sqrt(-beta_dens_der2(BETA_MODE,&used_data)/(2*PI)));
						//cout << marker << ": "<< prob << ";";

						gamma(marker) = dist.bernoulli(prob);   // Sample the inclusion based on the probability

						if(gamma(marker) == 1){
							new_xinit << BETA_MODE - 0.1 , BETA_MODE,  BETA_MODE+0.05, BETA_MODE + 0.1;
							assignArray(p_xinit,new_xinit);
							err = arms(xinit,ninit,&xl,&xr,beta_dens,&used_data,&convex,
									npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
							errorCheck(err);
							beta(marker) = xsamp[0];
							used_data.epsilon = used_data.epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new

						}else{
							beta(marker) = 0;
							//If beta is 0, then we don't need to do the residual update anymore
						}


					}

					// 3. Alpha
					xl = 0.0; xr = 400.0;
					new_xinit << (used_data.alpha)*0.5, used_data.alpha,  (used_data.alpha)*1.5, (used_data.alpha)*3;  // New values for abscissae evaluation
					assignArray(p_xinit,new_xinit);


					err = arms(xinit,ninit,&xl,&xr,alpha_dens,&used_data,&convex,
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
				left_trunc = data.left_trunc;


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

				/* Prior value selection for the variables */
				/* At the moment we set them to be uninformative */
				/* alpha */
				used_data.alpha_0 = alpha_0;
				used_data.kappa_0 = kappa_0;
				/* mu */
				used_data.sigma_mu = sigma_mu;
				/* sigma_b */
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



					/* 1. Mu */
					xl = -5; xr = 10.0;
					new_xinit << 0.95*mu, mu,  1.05*mu, 1.1*mu;  // New values for abscissae evaluation
					assignArray(p_xinit,new_xinit);
					used_data.epsilon = used_data.epsilon.array() + mu;//  we add the previous value
					used_data.epsilon_trunc = used_data.epsilon_trunc + mu;

					err = arms(xinit,ninit,&xl,&xr,mu_dens_ltrunc,&used_data,&convex,
							npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);

					errorCheck(err);
					mu = xsamp[0];
					used_data.epsilon = used_data.epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta
					used_data.epsilon_trunc = used_data.epsilon_trunc - mu;


					std::random_shuffle(markerI.begin(), markerI.end());

					// This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
					for (int j = 0; j < M; j++) {

						marker = markerI[j];

						// Using the preprocessed solution
						used_data.X_j = data.mappedZ.col(marker).cast<double>();

						used_data.epsilon = used_data.epsilon.array() + (used_data.X_j * beta(marker)).array();
						used_data.epsilon_trunc = used_data.epsilon_trunc.array() + (used_data.X_j * beta(marker)).array();


						/* Calculate the inclusion probability */
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
