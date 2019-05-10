#ifndef ANALYSISGRAPH_H
#define ANALYSISGRAPH_H

#include <cstddef>
#include <vector>

class BayesRBase;

class AnalysisGraph
{
public:
    AnalysisGraph(size_t maxParallel = 0);
    virtual ~AnalysisGraph();

    // Return true if the analysis calls BayesRBase::processColumnAsync
    virtual bool isAsynchronous() const = 0;

    virtual void exec(BayesRBase* bayes,
                      unsigned int numInds,
                      unsigned int numSnps,
                      const std::vector<unsigned int> &markerIndices) = 0;

protected:
    BayesRBase *m_bayes = nullptr;
    size_t m_maxParallel = 0; // Default to tbb::flow::unlimited
};

#endif // ANALYSISGRAPH_H
