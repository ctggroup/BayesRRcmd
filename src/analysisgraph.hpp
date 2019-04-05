#ifndef ANALYSISGRAPH_H
#define ANALYSISGRAPH_H

#include <cstddef>
#include <vector>

class BayesRBase;

class AnalysisGraph
{
public:
    AnalysisGraph(size_t maxParallel = 12);
    virtual ~AnalysisGraph();

    // Return true if the analysis calls BayesRBase::processColumnAsync
    virtual bool isAsynchronous() const = 0;

    virtual void exec(BayesRBase* bayes,
                      unsigned int numInds,
                      unsigned int numSnps,
                      const std::vector<unsigned int> &markerIndices) = 0;

protected:
    BayesRBase *m_bayes = nullptr;
    size_t m_maxParallel = 12;
};

#endif // ANALYSISGRAPH_H
