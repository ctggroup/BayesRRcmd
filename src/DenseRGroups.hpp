
#ifndef SRC_DENSERGROUPS_H_
#define SRC_DENSERGROUPS_H_

#include "BayesRGroupsBase.hpp"

class DenseRGroups : public BayesRGroupsBase
{

    friend class LimitSequenceGraph;
    friend class DenseParallelGraph;

public:
	explicit DenseRGroups(const Data *data, Options &opt);
    ~DenseRGroups() override;


    MarkerBuilder *markerBuilder() const override;

    void updateGlobal(Marker *marker, const double beta_old, const double beta,VectorXd& deltaEps) override;

	void updateMu(double old_mu,double N) override;

protected:
    void init(int K, unsigned int markerCount, unsigned int individualCount, unsigned int groupCount) override;

    void prepareForAnylsis() override;

    void readWithSharedLock(Marker *marker) override;
    
};

#endif /* SRC_DENSERGROUPS_H_ */

