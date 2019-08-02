#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "data.hpp"
#include "options.hpp"

#include <Eigen/Eigen>
#include <memory>

using namespace Eigen;

struct Kernel;
struct IndexEntry;
struct Marker;

class AnalysisGraph;
class MarkerBuilder;

struct AsyncResult {
    double betaOld = 0.0;
    double beta = 0.0;
    std::unique_ptr<VectorXd> deltaEpsilon;
};

class Analysis {
public:
    explicit Analysis(const Data *data, const Options &opt);
    virtual ~Analysis();

    virtual std::unique_ptr<Kernel> kernelForMarker(const Marker *marker) const = 0;

    virtual MarkerBuilder* markerBuilder() const = 0;
    virtual IndexEntry indexEntry(unsigned int i) const;
    virtual bool compressed() const;
    virtual unsigned char* compressedData() const;
    virtual std::string preprocessedFile() const;

    virtual int runGibbs(AnalysisGraph* analysis) = 0;

    virtual void processColumn(Kernel *kernel) = 0;
    virtual std::unique_ptr<AsyncResult> processColumnAsync(Kernel *kernel) = 0;

    virtual void updateGlobal(Kernel *kernel, const double beta_old, const double beta, const VectorXd& deltaEps) = 0;

protected:
    const Data *m_data = nullptr; // data matrices
    const Options &m_opt;
};

#endif // ANALYSIS_H
