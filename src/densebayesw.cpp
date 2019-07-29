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
#include "densebayesw.h"
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

DenseBayesW::DenseBayesW(Data &data, Options &opt, const long memPageSize)
: BayesWBase(data, opt, memPageSize)
{

}

namespace  {

/* Function to check if ARS resulted with error*/
inline void errorCheck(int err){
    if(err>0){
        cout << "Error code = " << err << endl;
        exit(1);
    }
}

struct beta_params {
    double alpha = 0;
    VectorXd epsilon;
    double sigma_b = 0;
    double sum_failure = 0;
    VectorXd X_j;

    VectorXd mixture_classes;
    int used_mixture = 0;
};

/* Function for the log density of beta: uses mixture component from the structure norm_data */
inline double beta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
    beta_params p = *(static_cast<beta_params *>(norm_data));

	y = -p.alpha * x * p.sum_failure - (((p.epsilon - p.X_j * x) * p.alpha).array() - EuMasc).exp().sum() -
			x * x / (2 * p.mixture_classes(p.used_mixture) * p.sigma_b) ;
	return y;
};

}


void DenseBayesW::sampleBeta(int marker)
{
    // Save the SNP effect column to a structure
    used_data.X_j = data.Z.col(marker).cast<double>();

    // Save sum(X_j*failure) to structure
    used_data_beta.sum_failure = sum_failure(marker);

    BayesWBase::sampleBeta(marker);
}

double DenseBayesW::calculateSumFailure(int marker)
{
    return ((data.Z.col(marker).cast<double>()).array() * used_data_alpha.failure_vector.array()).sum();
}

void DenseBayesW::preEstimateResidualUpdate(int marker)
{
    epsilon = epsilon.array() + (used_data.X_j * beta(marker)).array();
}

std::unique_ptr<gh_params> DenseBayesW::gaussHermiteParameters(int marker)
{
    return std::make_unique<dense_gh_params>(vi,used_data.X_j);
}

int DenseBayesW::estimateBeta(double *xinit, int ninit, double *xl, double *xr, double *convex, int npoint,
                              int dometrop, double *xprev, double *xsamp, int nsamp, double *qcent,
                              double *xcent, int ncent, int *neval)
{
    beta_params params;
    params.alpha = used_data.alpha;
    params.epsilon = epsilon;
    params.sigma_b = used_data_beta.sigma_b;
    params.sum_failure = used_data.sum_failure;
    params.X_j = used_data.X_j;
    params.mixture_classes = used_data.mixture_classes;
    params.used_mixture = used_data.used_mixture;

    return arms(xinit, ninit, xl, xr, beta_dens, &params, convex,
                npoint, dometrop, xprev, xsamp, nsamp, qcent, xcent, ncent, neval);
}

void DenseBayesW::postEstimateResidualUpdate(int marker)
{
    epsilon = epsilon - used_data.X_j * beta(marker); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new
}

double dense_gh_params::exponent_sum() const
{
    return (vi.array() * Xj.array() * Xj.array()).sum();
}

double dense_gh_params::integrand_adaptive(double s, double alpha, double dj, double sqrt_2Ck_sigmab) const
{
    //vi is a vector of exp(vi)
    double temp = -alpha *s*dj*sqrt_2Ck_sigmab + (vi.array()* (1 - (-Xj.array()*s*sqrt_2Ck_sigmab*alpha).exp() )).sum() -pow(s,2);
    return exp(temp);
}
