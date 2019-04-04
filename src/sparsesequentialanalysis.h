#ifndef SPARSESEQUENTIALANALYSIS_H
#define SPARSESEQUENTIALANALYSIS_H

#include "analysisgraph.hpp"

class SparseBayesRRG;

class SparseSequentialAnalysis : public AnalysisGraph
{
public:
    SparseSequentialAnalysis(SparseBayesRRG *bayes);

    void exec(unsigned int numKeptInds,
              unsigned int numIncdSnps,
              const std::vector<unsigned int> &markerIndices);

private:
    SparseBayesRRG *m_bayes;
};

#endif // SPARSESEQUENTIALANALYSIS_H
