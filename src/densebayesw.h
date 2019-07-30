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

struct DenseGaussMarker : public GaussMarker {
    DenseGaussMarker(const VectorXd &vi, const VectorXd &Zj, int i)
        : GaussMarker(i), vi(vi), Zj(Zj) {}

    const VectorXd &vi;
    const VectorXd &Zj;

    double exponent_sum() const override;
    double integrand_adaptive(double s, double alpha, double sqrt_2Ck_sigmab) const override;
};

class DenseBayesW : public BayesWBase
{
protected:
    VectorXd Z_j;

public:
    DenseBayesW(Data &data, Options &opt, const long memPageSize);

protected:
    void sampleBeta(int marker) override;

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


#endif /* DENSEBAYESW_HPP_ */
