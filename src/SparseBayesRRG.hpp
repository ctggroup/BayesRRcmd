#ifndef SRC_SPARSEBAYESRRG_H_
#define SRC_SPARSEBAYESRRG_H_

#include "BayesRBase.hpp"

#include <shared_mutex>

class EigenSparseData;

class SparseBayesRRG : public BayesRBase
{
    friend class SparseParallelGraph;

public:
    explicit SparseBayesRRG(const EigenSparseData *m_data, Options &m_opt);
    ~SparseBayesRRG() override;

    void processColumn(unsigned int marker);
    std::tuple<double, double> processColumnAsync(unsigned int marker);
    void updateGlobal(const unsigned int marker, double beta_old, double beta);

protected:
    const EigenSparseData *m_sparseData;

    VectorXd m_asyncEpsilon;

    double m_asyncEpsilonSum = 0.0;

    mutable std::shared_mutex m_mutex;
    mutable std::mutex m_rngMutex;

    void init(int K, unsigned int markerCount, unsigned int individualCount) override;
    void prepareForAnylsis() override;
};

#endif /* SRC_SPARSEBAYESRRG_H_ */
