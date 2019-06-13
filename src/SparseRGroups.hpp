#ifndef SRC_SPARSERGROUPS_H_
#define SRC_SPARSERGROUPS_H_

#include "BayesRGroupsBase.hpp"
#include <Eigen/Eigen>
class Data;

class SparseRGroups : public BayesRGroupsBase
{
    friend class SparseParallelGraph;

public:
    explicit SparseRGroups(const Data *m_data, Options &m_opt);
    ~SparseRGroups() override;

    MarkerBuilder *markerBuilder() const override;

   void updateGlobal(Marker *marker, const double beta_old, const double beta,VectorXd& deltaEps ) override;
   void updateMu(double old_mu,double N);

protected:
    double m_asyncEpsilonSum = 0.0;

    VectorXd m_ones;

    void init(int K, unsigned int markerCount, unsigned int individualCount, unsigned int groupCount) override;
    void prepareForAnylsis() override;

    void prepare(Marker *marker) override;
    void readWithSharedLock(Marker *marker) override;
    void writeWithUniqueLock(Marker *marker) override;
   
};

#endif /* SRC_SPARSERGROUPS_H_ */
