/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef SPARSEBAYESW_HPP_
#define SPARSEBAYESW_HPP_

#include "bayeswbase.h"
#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

struct sparse_gh_params : public gh_params {
    sparse_gh_params(double vi_sum, double vi_2, double vi_1, double vi_0, double mean, double sd, double mean_sd_ratio) :
        vi_sum(vi_sum), vi_2(vi_2), vi_1(vi_1), vi_0(vi_0), mean(mean), sd(sd), mean_sd_ratio(mean_sd_ratio) {}

    double vi_sum = 0;
    double vi_2 = 0;
    double vi_1 = 0;
    double vi_0 = 0;
    double mean = 0;
    double sd = 0;
    double mean_sd_ratio = 0;

    double exponent_sum() const override;
    double integrand_adaptive(double s ,double alpha, double dj, double sqrt_2Ck_sigmab) const override;
};

class SparseBayesW : public BayesWBase
{
public:
    SparseBayesW(Data &data, Options &opt, const long memPageSize);

protected:
    void sampleBeta(int marker) override;

    double calculateSumFailure(int marker) override;

    void preEstimateResidualUpdate(int marker) override;

    std::unique_ptr<gh_params> gaussHermiteParameters(int marker) override;

    int estimateBeta(int marker, double *xinit, int ninit, double *xl, double *xr, const beta_params params,
                      double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                      int nsamp, double *qcent, double *xcent,
                      int ncent, int *neval) override;

    void postEstimateResidualUpdate(int marker) override;
};


#endif /* SPARSEBAYESW_HPP_ */
