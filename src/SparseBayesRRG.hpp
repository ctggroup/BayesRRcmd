#ifndef SRC_SPARSEBAYESRRG_H_
#define SRC_SPARSEBAYESRRG_H_

#include "BayesRBase.hpp"

class SparseData;

class SparseBayesRRG : public BayesRBase
{
public:
    explicit SparseBayesRRG(const SparseData *m_data, Options &m_opt);
    ~SparseBayesRRG() override;

    void processColumn(unsigned int marker);

protected:
    const SparseData *m_sparseData;

    double m_epsilonSum = 0.0;
    VectorXd m_ones;

    // Helper references into our sparse data
    const VectorXd &m_means;
    const VectorXd &m_sds;
    const VectorXd &m_sqrdZ;
    const VectorXd &m_Zsum;

    void init(int K, unsigned int markerCount, unsigned int individualCount) override;

    double dot(const unsigned int marker, const VectorXd &epsilon, const double sd) const;
};

#endif /* SRC_SPARSEBAYESRRG_H_ */
