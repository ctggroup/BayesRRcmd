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

    void processColumn(unsigned int marker, const Map<VectorXd> &Cx);
    std::tuple<double, double> processColumnAsync(unsigned int marker, const Map<VectorXd> &Cx);

    void updateGlobal(double beta_old, double beta, const Map<VectorXd> &Cx);
    void updateGlobal(Marker *marker, const double beta_old, const double beta) override;

protected:
    void init(int K, unsigned int markerCount, unsigned int individualCount) override;
    void prepareForAnylsis() override;

    void readWithSharedLock(Marker *marker) override;
};

#endif /* SRC_DENSEBAYESRRMZ_H_ */
