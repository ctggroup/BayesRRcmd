/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef DENSEBAYESW_HPP_
#define DENSEBAYESW_HPP_

#include "bayeswbase.h"

class DenseBayesW : public BayesWBase
{
public:
    DenseBayesW(const Data *data, const Options *opt, const long m_memPageSize);

    std::unique_ptr<Kernel> kernelForMarker(const ConstMarkerPtr &marker) const override;
    MarkerBuilder *markerBuilder() const override;

protected:
    int estimateBeta(const BayesWKernel *kernel, const std::shared_ptr<VectorXd> &epsilon, double *xinit, int ninit, double *xl, double *xr, const beta_params params,
                      double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
                      int nsamp, double *qcent, double *xcent,
                      int ncent, int *neval) override;
};


#endif /* DENSEBAYESW_HPP_ */
