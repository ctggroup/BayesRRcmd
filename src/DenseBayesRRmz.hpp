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
    explicit DenseBayesRRmz(const Data *m_data, Options &m_opt);
    ~DenseBayesRRmz() override;

    MarkerBuilder *markerBuilder() const override;

  void updateGlobal(Marker *marker, const double beta_old, const double beta,VectorXd& deltaEps) override;

protected:
    void init(int K, unsigned int markerCount, unsigned int individualCount) override;
    void prepareForAnylsis() override;

    void readWithSharedLock(Marker *marker) override;
};

#endif /* SRC_DENSEBAYESRRMZ_H_ */
