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

struct SparseGaussMarker : public GaussMarker {
    SparseGaussMarker(int i) : GaussMarker(i) {}

    double vi_sum = 0;
    double vi_2 = 0;
    double vi_1 = 0;
    double vi_0 = 0;
    double mean = 0;
    double sd = 0;
    double mean_sd_ratio = 0;

    double exponent_sum() const override;
    double integrand_adaptive(double s , double alpha, double sqrt_2Ck_sigmab) const override;
};

class SparseBayesW : public BayesWBase
{
public:
    SparseBayesW(Data &data, Options &opt, const long memPageSize);

protected:
    double calculateSumFailure(int marker) override;

    std::unique_ptr<GaussMarker> buildMarker(int i) override;
    void prepare(GaussMarker *marker) override;

    void preEstimateResidualUpdate(const GaussMarker *marker) override;

    int estimateBeta(const GaussMarker *marker, double *xinit, int ninit, double *xl, double *xr, const beta_params params,
                      double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                      int nsamp, double *qcent, double *xcent,
                      int ncent, int *neval) override;

    void postEstimateResidualUpdate(const GaussMarker *marker) override;
};


#endif /* SPARSEBAYESW_HPP_ */
