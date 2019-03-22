#ifndef SRC_SPARSEBAYESRRG_H_
#define SRC_SPARSEBAYESRRG_H_

#include "BayesRBase.hpp"

#include <shared_mutex>

class SparseData;

class SparseBayesRRG : public BayesRBase
{
    friend class SparseParallelGraph;

public:
    explicit SparseBayesRRG(const SparseData *m_data, Options &m_opt);
    ~SparseBayesRRG() override;

    void processColumn(unsigned int marker);
    std::tuple<double, double> processColumnAsync(unsigned int marker);
    void updateGlobal(const unsigned int marker, double beta_old, double beta);

protected:
    const SparseData *m_sparseData;

    VectorXd m_asyncEpsilon;

    double m_asyncEpsilonSum = 0.0;

    VectorXd m_ones;

    // Helper references into our sparse data
    const VectorXd &m_means;
    const VectorXd &m_sds;
    const VectorXd &m_sqrdZ;
    const VectorXd &m_Zsum;

    mutable std::shared_mutex m_mutex;
    mutable std::mutex m_rngMutex;

    void init(int K, unsigned int markerCount, unsigned int individualCount) override;
    void prepareForAnylsis() override;

    double computeNum(const unsigned int marker, const double beta_old, const VectorXd &epsilon) const;
    double dot(const unsigned int marker, const VectorXd &epsilon, const double sd) const;
    VectorXd computeEpsilonUpdate(const unsigned int marker, const double beta_old, const double beta) const;
    double computeEpsilonSumUpdate(const unsigned int marker, const double beta_old, const double beta) const;
};

#endif /* SRC_SPARSEBAYESRRG_H_ */
