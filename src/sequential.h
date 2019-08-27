#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include "analysisgraph.hpp"

class Sequential : public AnalysisGraph
{
public:
    bool isAsynchronous() const override { return false; }
    void exec(Analysis *analysis,
              unsigned int numInds,
              unsigned int numSnps,
              const std::vector<unsigned int> &markerIndices) override;
};

#endif
