/*
 * BayesRRm.h
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#ifndef SRC_DENSEBAYESRRMZ_H_
#define SRC_DENSEBAYESRRMZ_H_

#include "BayesRBase.hpp"

#include <shared_mutex>

class DenseBayesRRmz : public BayesRBase
{
    friend class LimitSequenceGraph;
    friend class ParallelGraph;

public:
    explicit DenseBayesRRmz(const Data *m_data, Options &m_opt);
    ~DenseBayesRRmz() override;

    void processColumn(unsigned int marker, const Map<VectorXd> &Cx);
    std::tuple<double, double> processColumnAsync(unsigned int marker, const Map<VectorXd> &Cx);
    void updateGlobal(double beta_old, double beta, const Map<VectorXd> &Cx);

protected:
    mutable std::shared_mutex m_mutex;
    mutable std::mutex m_rngMutex;
};

#endif /* SRC_DENSEBAYESRRMZ_H_ */
