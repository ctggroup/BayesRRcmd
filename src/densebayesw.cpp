/*
 * BayesW.cpp
 *
 *  Created on: 26 Nov 2018
 *  Author: Sven Erik Ojavee
 *  Last changes: 22 Feb 2019
 */

#include "densebayesw.h"
#include "densebayeswkernel.h"
#include "BayesW_arms.h"

/* Pre-calculate used constants */
#define EuMasc 0.577215664901532

namespace  {

struct dense_beta_params : public beta_params {
    dense_beta_params(const beta_params &params) : beta_params(params) {}
    std::shared_ptr<VectorXd> epsilon = nullptr;
    std::shared_ptr<Map<VectorXd>> Cx = nullptr;
};

/* Function for the log density of beta: uses mixture component from the structure norm_data */
inline double beta_dens(double x, void *norm_data)
/* We are sampling beta (denoted by x here) */
{
    /* In C++ we need to do a static cast for the void data */
    dense_beta_params p = *(static_cast<dense_beta_params *>(norm_data));

    return -p.alpha * x * p.sum_failure - (((*p.epsilon - *p.Cx * x) * p.alpha).array() - EuMasc).exp().sum() -
            x * x / (2 * p.used_mixture * p.sigma_b) ;
};

}

DenseBayesW::DenseBayesW(const Data *data, const Options *opt, const long memPageSize)
: BayesWBase(data, opt, memPageSize)
{

}

std::unique_ptr<Kernel> DenseBayesW::kernelForMarker(const ConstMarkerPtr &marker) const
{
    const auto denseMarker = dynamic_pointer_cast<const DenseMarker>(marker);
    assert(denseMarker);
    return std::make_unique<DenseBayesWKernel>(denseMarker);
}

MarkerBuilder *DenseBayesW::markerBuilder() const
{
    return builderForType(PreprocessDataType::Dense);
}

int DenseBayesW::estimateBeta(const BayesWKernel *kernel, const std::shared_ptr<VectorXd> &epsilon, double *xinit, int ninit, double *xl, double *xr, const beta_params params, double *convex, int npoint,
                              int dometrop, double *xprev, double *xsamp, int nsamp, double *qcent,
                              double *xcent, int ncent, int *neval)
{
    const auto* denseMarker = dynamic_cast<const DenseMarker*>(kernel->marker.get());
    assert(denseMarker);

    dense_beta_params dense_params {params};
    dense_params.epsilon = epsilon;
    dense_params.Cx = denseMarker->Cx;

    return arms(xinit, ninit, xl, xr, beta_dens, &dense_params, convex,
                npoint, dometrop, xprev, xsamp, nsamp, qcent, xcent, ncent, neval);
}
