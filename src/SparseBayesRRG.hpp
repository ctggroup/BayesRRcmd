#ifndef SRC_SPARSEBAYESRRG_H_
#define SRC_SPARSEBAYESRRG_H_

#include "BayesRBase.hpp"

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
    void updateGlobal(Marker *marker, const double beta_old, const double beta) override;

protected:
    const SparseData *m_sparseData;

    double m_asyncEpsilonSum = 0.0;

    VectorXd m_ones;

    void init(int K, unsigned int markerCount, unsigned int individualCount) override;
    void prepareForAnylsis() override;

    void prepare(Marker *marker) override;
    void readWithSharedLock(Marker *marker) override;
    void writeWithUniqueLock(Marker *marker) override;
};

#endif /* SRC_SPARSEBAYESRRG_H_ */
