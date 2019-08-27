/*
 * BayesRRm.h
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#ifndef SRC_DENSEBAYESRRMZ_H_
#define SRC_DENSEBAYESRRMZ_H_

#include "BayesRBase.hpp"

class DenseBayesRRmz : public BayesRBase
{
    friend class LimitSequenceGraph;
    friend class DenseParallelGraph;

public:
    explicit DenseBayesRRmz(const Data *data, const Options *opt);
    ~DenseBayesRRmz() override;

    std::unique_ptr<Kernel> kernelForMarker(const ConstMarkerPtr &marker) const override;
    MarkerBuilder *markerBuilder() const override;

    void updateMu(double old_mu,double N) override;
};

#endif /* SRC_DENSEBAYESRRMZ_H_ */
