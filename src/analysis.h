#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "asyncresult.h"
#include "common.h"
#include "data.hpp"
#include "options.hpp"

#include <Eigen/Eigen>
#include <memory>

using namespace Eigen;

class AnalysisGraph;
class MarkerBuilder;

class Analysis {
public:
    explicit Analysis(const Data *data, const Options *opt);
    virtual ~Analysis();

    virtual std::unique_ptr<Kernel> kernelForMarker(const ConstMarkerPtr &marker) const = 0;

    virtual MarkerBuilder* markerBuilder() const = 0;
    virtual IndexEntry indexEntry(unsigned int i) const;
    virtual bool compressed() const;
    virtual unsigned char* compressedData() const;
    virtual std::string preprocessedFile() const;

    virtual int runGibbs(AnalysisGraph* analysis);
    virtual int runGibbs(AnalysisGraph* analysis, const std::vector<unsigned int> &markers);
    virtual int runGibbs(AnalysisGraph* analysis, std::vector<unsigned int> &&markers) = 0;

    // LimitSeqeunceGraph
    virtual void processColumn(const KernelPtr &kernel) = 0;

    // ParallelGraph
    virtual std::unique_ptr<AsyncResult> processColumnAsync(const KernelPtr &kernel) = 0;
    virtual void doThreadSafeUpdates(const ConstAsyncResultPtr& result) = 0;
    virtual void updateGlobal(const KernelPtr& kernel,
                              const ConstAsyncResultPtr& result) = 0;

protected:
    const Data *m_data = nullptr; // data matrices
    const Options *m_opt;

    // MPI
    int m_rank = 0;
    int m_worldSize = 1;
};

#endif // ANALYSIS_H
