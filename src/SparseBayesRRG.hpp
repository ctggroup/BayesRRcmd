#ifndef SRC_SPARSEBAYESRRG_H_
#define SRC_SPARSEBAYESRRG_H_

#include "BayesRBase.hpp"
#include <Eigen/Eigen>
class Data;

class SparseBayesRRG : public BayesRBase
{
    friend class SparseParallelGraph;

public:
    explicit SparseBayesRRG(const Data *data, const Options *opt);
    ~SparseBayesRRG() override;

    std::unique_ptr<Kernel> kernelForMarker(const ConstMarkerPtr &marker) const override;
    MarkerBuilder *markerBuilder() const override;

    void updateGlobal(const KernelPtr& kernel, const ConstAsyncResultPtr &result) override;
    void updateMu(double old_mu,double N);
protected:
    VectorXd m_ones;

    void init(int K, unsigned int markerCount, unsigned int individualCount) override;

    void prepare(BayesRKernel *kernel) override;
    void readWithSharedLock(BayesRKernel *kernel) override;
    void writeWithUniqueLock(BayesRKernel *kernel) override;
   
};

#endif /* SRC_SPARSEBAYESRRG_H_ */
