/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef DENSEBAYESW_HPP_
#define DENSEBAYESW_HPP_

#include "bayeswbase.h"
#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

struct dense_gh_params : public gh_params {
    dense_gh_params(const VectorXd &vi, const VectorXd &Xj) : vi(vi), Xj(Xj) {}
    const VectorXd &vi;
    const VectorXd &Xj;

    double exponent_sum() const override;
    double integrand_adaptive(double s, double alpha, double dj, double sqrt_2Ck_sigmab) const override;
};

class DenseBayesW : public BayesWBase
{
public:
    DenseBayesW(Data &data, Options &opt, const long memPageSize);

protected:
    void sampleBeta(int marker) override;

    double calculateSumFailure(int marker) override;

    void preEstimateResidualUpdate(int marker) override;

    std::unique_ptr<gh_params> gaussHermiteParameters(int marker) override;

    int estimateBeta (double *xinit, int ninit, double *xl, double *xr,
                      double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                      int nsamp, double *qcent, double *xcent,
                      int ncent, int *neval) override;

    void postEstimateResidualUpdate(int marker) override;

};


#endif /* DENSEBAYESW_HPP_ */
