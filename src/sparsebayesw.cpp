/*
 * BayesW.cpp
 *
 *  Created on: 26 Nov 2018
 *  Author: Sven Erik Ojavee
 *  Last changes: 22 Feb 2019
 */

#include "raggedbayeswkernel.h"
#include "sparsebayesw.h"
#include "BayesW_arms.h"

namespace {

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

SparseBayesW::SparseBayesW(const Data *data, const Options *opt, const long memPageSize)
: BayesWBase (data, opt, memPageSize)
{

}

std::unique_ptr<Kernel> SparseBayesW::kernelForMarker(const ConstMarkerPtr &marker) const
{
    switch (m_opt->preprocessDataType) {
    case PreprocessDataType::SparseRagged:
    {
        const auto raggedSparseMarker = dynamic_pointer_cast<const RaggedSparseMarker>(marker);
        assert(raggedSparseMarker);
        return std::make_unique<RaggedBayesWKernel>(raggedSparseMarker);
    }

    default:
        std::cerr << "SparseBayesW::kernelForMarker - unsupported type: "
                  << m_opt->preprocessDataType
                  << std::endl;
    }

    return {};
}

MarkerBuilder *SparseBayesW::markerBuilder() const
{
    switch (m_opt->preprocessDataType) {
    case PreprocessDataType::SparseRagged:
        return builderForType(PreprocessDataType::SparseRagged);

    default:
        std::cerr << "SparseBayesW::markerBuilder - unsupported type: "
                  << m_opt->preprocessDataType
                  << std::endl;
    }

    return nullptr;
}

int SparseBayesW::estimateBeta(const BayesWKernel *kernel, const std::shared_ptr<VectorXd> &epsilon, double *xinit, int ninit, double *xl, double *xr, const beta_params params, double *convex, int npoint,
                               int dometrop, double *xprev, double *xsamp, int nsamp, double *qcent,
                               double *xcent, int ncent, int *neval)
 {
    (void) epsilon; // Unused

    const auto* raggedKernel = dynamic_cast<const RaggedBayesWKernel*>(kernel);
    assert(raggedKernel);

    const auto* raggedMarker = dynamic_cast<const RaggedSparseMarker*>(kernel->marker.get());
    assert(raggedMarker);

    sparse_beta_params sparse_params {params};
    sparse_params.mean_sd_ratio = raggedMarker->mean / raggedMarker->sd;
    sparse_params.sd = raggedMarker->sd;
    sparse_params.vi_0 = raggedKernel->vi_0;
    sparse_params.vi_1 = raggedKernel->vi_1;
    sparse_params.vi_2 = raggedKernel->vi_2;

     return arms(xinit, ninit, xl, xr, beta_dens, &sparse_params, convex,
                 npoint, dometrop, xprev, xsamp, nsamp, qcent, xcent, ncent, neval);
}
