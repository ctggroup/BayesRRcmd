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

struct dense_beta_params : public beta_params {
    dense_beta_params(const beta_params &params) : beta_params(params) {}
    VectorXd epsilon;
    VectorXd Z_j;
};

/* Function for the log density of beta: uses mixture component from the structure norm_data */
inline double beta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
	double y;

	/* In C++ we need to do a static cast for the void data */
    dense_beta_params p = *(static_cast<dense_beta_params *>(norm_data));

    y = -p.alpha * x * p.sum_failure - (((p.epsilon - p.Z_j * x) * p.alpha).array() - EuMasc).exp().sum() -
            x * x / (2 * p.used_mixture * p.sigma_b) ;
	return y;
};

}


void DenseBayesW::sampleBeta(int marker)
{
    // Save the SNP effect column to a structure
    Z_j = data.Z.col(marker).cast<double>();

    BayesWBase::sampleBeta(marker);
}

double DenseBayesW::calculateSumFailure(int marker)
{
    return ((data.Z.col(marker).cast<double>()).array() * failure_vector.array()).sum();
}

std::unique_ptr<GaussMarker> DenseBayesW::buildMarker(int i)
{
    return std::make_unique<DenseGaussMarker>(vi, Z_j, i);
}

void DenseBayesW::prepare(GaussMarker *marker)
{
    BayesWBase::prepare(marker);

    if (auto* denseMarker = dynamic_cast<DenseGaussMarker*>(marker)) {

    }
}

void DenseBayesW::preEstimateResidualUpdate(const GaussMarker *marker)
{
    const auto* denseMarker = dynamic_cast<const DenseGaussMarker*>(marker);
    assert(denseMarker);

    epsilon = epsilon.array() + (Z_j * beta(denseMarker->i)).array();
}

int DenseBayesW::estimateBeta(const GaussMarker *marker, double *xinit, int ninit, double *xl, double *xr, const beta_params params, double *convex, int npoint,
                              int dometrop, double *xprev, double *xsamp, int nsamp, double *qcent,
                              double *xcent, int ncent, int *neval)
{
    (void) marker; // Unused

    dense_beta_params dense_params {params};
    dense_params.epsilon = epsilon;
    dense_params.sum_failure = marker->sum_failure;
    dense_params.Z_j = Z_j;

    return arms(xinit, ninit, xl, xr, beta_dens, &dense_params, convex,
                npoint, dometrop, xprev, xsamp, nsamp, qcent, xcent, ncent, neval);
}

void DenseBayesW::postEstimateResidualUpdate(const GaussMarker *marker)
{
    const auto* denseMarker = dynamic_cast<const DenseGaussMarker*>(marker);
    assert(denseMarker);

    epsilon = epsilon - Z_j * beta(denseMarker->i); //now epsilon contains Y-mu - X*beta+ X.col(marker)*beta(marker)_old- X.col(marker)*beta(marker)_new
}

double DenseGaussMarker::exponent_sum() const
{
    return (vi.array() * Zj.array() * Zj.array()).sum();
}

double DenseGaussMarker::integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const
{
    //vi is a vector of exp(vi)
    double temp = -alpha *s*sum_failure*sqrt_2Ck_sigmab + (vi.array()* (1 - (-Zj.array()*s*sqrt_2Ck_sigmab*alpha).exp() )).sum() -pow(s,2);
    return exp(temp);
}
