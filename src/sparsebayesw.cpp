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

struct sparse_beta_params : beta_params {
    sparse_beta_params(const beta_params &params) : beta_params(params) {}
    double mean_sd_ratio = 0;
    double sd = 0;

    double vi_0 = 0;
    double vi_1 = 0;
    double vi_2 = 0;
};

/* Sparse version for function for the log density of beta: uses mixture component from the structure norm_data */
inline double beta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;
	/* In C++ we need to do a static cast for the void data */
    sparse_beta_params p = *(static_cast<sparse_beta_params *>(norm_data));

	y = -p.alpha * x * p.sum_failure -
			exp(p.alpha*x*p.mean_sd_ratio)* (p.vi_0 + p.vi_1 * exp(-p.alpha*x/p.sd) + p.vi_2 * exp(-2*p.alpha*x/p.sd))
            -x * x / (2 * p.used_mixture * p.sigma_b) ;
	return y;
};

}

double SparseBayesW::calculateSumFailure(int marker)
{
    std::vector<int> oneIndices = data.Zones[marker]; //Take the vector of indices
    std::vector<int> twoIndices = data.Ztwos[marker]; //Take the vector of indices

    int temp_sum = 0;
    for(int i=0; i < oneIndices.size(); i++){
        temp_sum += failure_vector(oneIndices[i]);
    }
    for(int i=0; i < twoIndices.size(); i++){
        temp_sum += 2*failure_vector(twoIndices[i]);
    }

    return (temp_sum - data.means(marker) * failure_vector.array().sum()) / data.sds(marker);
}

std::unique_ptr<GaussMarker> SparseBayesW::buildMarker(int i)
{
    auto marker = std::make_unique<SparseGaussMarker>(i);

    marker->vi_sum = vi.sum();
    marker->vi_2 = vi(data.Ztwos[i]).sum();
    marker->vi_1 = vi(data.Zones[i]).sum();
    marker->vi_0 = marker->vi_sum - marker->vi_1 - marker->vi_2;
    marker->mean = data.means(i);
    marker->sd = data.sds(i);
    marker->mean_sd_ratio = data.mean_sd_ratio(i);

    return std::move(marker);
}

void SparseBayesW::prepare(GaussMarker *marker)
{
    BayesWBase::prepare(marker);

    if (auto* sparseMarker = dynamic_cast<SparseGaussMarker*>(marker)) {

    }
}

void SparseBayesW::preEstimateResidualUpdate(const GaussMarker *marker)
{
    const auto* sparseMarker = dynamic_cast<const SparseGaussMarker*>(marker);
    assert(sparseMarker);

    epsilon = epsilon.array() - sparseMarker->mean_sd_ratio * beta(sparseMarker->i);  //Adjust for every memeber
    //And adjust even further for specific 1 and 2 allele values
    for(int i=0; i < data.Zones[sparseMarker->i].size(); i++){
        epsilon[data.Zones[sparseMarker->i][i]] += beta(sparseMarker->i)/sparseMarker->sd;
    }
    for(int i=0; i < data.Ztwos[sparseMarker->i].size(); i++){
        epsilon[data.Ztwos[sparseMarker->i][i]] += 2*beta(sparseMarker->i)/sparseMarker->sd;
    }
}

int SparseBayesW::estimateBeta(const GaussMarker *marker, double *xinit, int ninit, double *xl, double *xr, const beta_params params, double *convex, int npoint,
                               int dometrop, double *xprev, double *xsamp, int nsamp, double *qcent,
                               double *xcent, int ncent, int *neval)
 {
    const auto* sparseMarker = dynamic_cast<const SparseGaussMarker*>(marker);
    assert(sparseMarker);

    sparse_beta_params sparse_params {params};
    sparse_params.mean_sd_ratio = sparseMarker->mean_sd_ratio;
    sparse_params.sd = sparseMarker->sd;
    sparse_params.sum_failure = sparseMarker->sum_failure;
    sparse_params.vi_0 = sparseMarker->vi_0;
    sparse_params.vi_1 = sparseMarker->vi_1;
    sparse_params.vi_2 = sparseMarker->vi_2;

     return arms(xinit, ninit, xl, xr, beta_dens, &sparse_params, convex,
                 npoint, dometrop, xprev, xsamp, nsamp, qcent, xcent, ncent, neval);
}

void SparseBayesW::postEstimateResidualUpdate(const GaussMarker *marker)
{
    const auto* sparseMarker = dynamic_cast<const SparseGaussMarker*>(marker);
    assert(sparseMarker);

    epsilon = epsilon.array() + sparseMarker->mean_sd_ratio * beta(sparseMarker->i);  //Adjust for every memeber
    //And adjust even further for specific 1 and 2 allele values
    for(int i=0; i < data.Zones[sparseMarker->i].size(); i++){
        epsilon[data.Zones[sparseMarker->i][i]] -= beta(sparseMarker->i)/sparseMarker->sd;
    }
    for(int i=0; i < data.Ztwos[sparseMarker->i].size(); i++){
        epsilon[data.Ztwos[sparseMarker->i][i]] -= 2*beta(sparseMarker->i)/sparseMarker->sd;
    }
}

double SparseGaussMarker::exponent_sum() const
{
    return (vi_1 * (1 - 2 * mean) + 4 * (1-mean) * vi_2 + vi_sum * mean * mean) /(sd*sd);
}

double SparseGaussMarker::integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const
{
    //vi is a vector of exp(vi)
    double temp = -alpha *s*sum_failure*sqrt_2Ck_sigmab +
            vi_sum - exp(alpha*mean_sd_ratio*s*sqrt_2Ck_sigmab) *
            (vi_0 + vi_1 * exp(-alpha * s*sqrt_2Ck_sigmab/sd) + vi_2* exp(-2 * alpha * s*sqrt_2Ck_sigmab/sd))
            -pow(s,2);
    return exp(temp);
}
