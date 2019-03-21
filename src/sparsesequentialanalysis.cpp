#include "sparsesequentialanalysis.h"

#include "SparseBayesRRG.hpp"

SparseSequentialAnalysis::SparseSequentialAnalysis(SparseBayesRRG *bayes)
    : AnalysisGraph(1)
    , m_bayes(bayes)
{

}

void SparseSequentialAnalysis::exec(unsigned int numKeptInds,
                                    unsigned int numIncdSnps,
                                    const std::vector<unsigned int> &markerIndices)
{
    (void)numKeptInds; // Unused

    for (unsigned int i = 0; i < numIncdSnps; ++i) {
        m_bayes->processColumn(markerIndices[i]);
    }
}
