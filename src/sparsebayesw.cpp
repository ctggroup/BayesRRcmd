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
#include "sparsebayesw.h"
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

SparseBayesW::SparseBayesW(Data &data, Options &opt, const long memPageSize)
: BayesWBase (data, opt, memPageSize)
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

/* Sparse version for function for the log density of beta: uses mixture component from the structure norm_data */
inline double beta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;
	/* In C++ we need to do a static cast for the void data */
	pars_beta_sparse p = *(static_cast<pars_beta_sparse *>(norm_data));

	y = -p.alpha * x * p.sum_failure -
			exp(p.alpha*x*p.mean_sd_ratio)* (p.vi_0 + p.vi_1 * exp(-p.alpha*x/p.sd) + p.vi_2 * exp(-2*p.alpha*x/p.sd))
			-x * x / (2 * p.mixture_classes(p.used_mixture) * p.sigma_b) ;
	return y;
};

}

void SparseBayesW::sampleBeta(int marker)
{
    //Save sum(X_j*failure) to structure
    used_data_beta.sum_failure = sum_failure(marker);

    used_data_beta.mean = data.means(marker);
    used_data_beta.sd = data.sds(marker);
    used_data_beta.mean_sd_ratio = data.mean_sd_ratio(marker);

    BayesWBase::sampleBeta(marker);
}

double SparseBayesW::calculateSumFailure(int marker)
{
    std::vector<int> oneIndices = data.Zones[marker]; //Take the vector of indices
    std::vector<int> twoIndices = data.Ztwos[marker]; //Take the vector of indices

    int temp_sum = 0;
    for(int i=0; i < oneIndices.size(); i++){
        temp_sum += used_data_alpha.failure_vector(oneIndices[i]);
    }
    for(int i=0; i < twoIndices.size(); i++){
        temp_sum += 2*used_data_alpha.failure_vector(twoIndices[i]);
    }

    return (temp_sum - data.means(marker) * used_data_alpha.failure_vector.array().sum()) / data.sds(marker);
}

void SparseBayesW::preEstimateResidualUpdate(int marker)
{
    epsilon = epsilon.array() - used_data_beta.mean_sd_ratio * beta(marker);  //Adjust for every memeber
    //And adjust even further for specific 1 and 2 allele values
    for(int i=0; i < data.Zones[marker].size(); i++){
        epsilon[data.Zones[marker][i]] += beta(marker)/used_data_beta.sd;
    }
    for(int i=0; i < data.Ztwos[marker].size(); i++){
        epsilon[data.Ztwos[marker][i]] += 2*beta(marker)/used_data_beta.sd;
    }
}

std::unique_ptr<gh_params> SparseBayesW::gaussHermiteParameters(int marker)
{
    double vi_sum = vi.sum();
    used_data_beta.vi_2 = vi(data.Ztwos[marker]).sum();
    used_data_beta.vi_1 = vi(data.Zones[marker]).sum();
    used_data_beta.vi_0 = vi_sum - used_data_beta.vi_1 - used_data_beta.vi_2;

    return std::make_unique<sparse_gh_params>(vi_sum, used_data_beta.vi_2, used_data_beta.vi_1,
                                              used_data_beta.vi_0, data.means(marker),data.sds(marker),data.mean_sd_ratio(marker));
}

int SparseBayesW::estimateBeta(double *xinit, int ninit, double *xl, double *xr, double *convex, int npoint,
                               int dometrop, double *xprev, double *xsamp, int nsamp, double *qcent,
                               double *xcent, int ncent, int *neval)
 {
     return arms(xinit, ninit, xl, xr, beta_dens, &used_data_beta, convex,
                 npoint, dometrop, xprev, xsamp, nsamp, qcent, xcent, ncent, neval);
}

void SparseBayesW::postEstimateResidualUpdate(int marker)
{
    epsilon = epsilon.array() + used_data_beta.mean_sd_ratio * beta(marker);  //Adjust for every memeber
    //And adjust even further for specific 1 and 2 allele values
    for(int i=0; i < data.Zones[marker].size(); i++){
        epsilon[data.Zones[marker][i]] -= beta(marker)/used_data_beta.sd;
    }
    for(int i=0; i < data.Ztwos[marker].size(); i++){
        epsilon[data.Ztwos[marker][i]] -= 2*beta(marker)/used_data_beta.sd;
    }
}

double sparse_gh_params::exponent_sum() const
{
    return (vi_1 * (1 - 2 * mean) + 4 * (1-mean) * vi_2 + vi_sum * mean * mean) /(sd*sd);
}

double sparse_gh_params::integrand_adaptive(double s, double alpha, double dj, double sqrt_2Ck_sigmab) const
{
    //vi is a vector of exp(vi)
    double temp = -alpha *s*dj*sqrt_2Ck_sigmab +
            vi_sum - exp(alpha*mean_sd_ratio*s*sqrt_2Ck_sigmab) *
            (vi_0 + vi_1 * exp(-alpha * s*sqrt_2Ck_sigmab/sd) + vi_2* exp(-2 * alpha * s*sqrt_2Ck_sigmab/sd))
            -pow(s,2);
    return exp(temp);
}
