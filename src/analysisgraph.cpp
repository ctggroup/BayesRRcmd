#include "analysisgraph.h"

AnalysisGraph::AnalysisGraph(BayesRRmz *bayes, size_t maxParallel)
    : m_bayes(bayes)
    , m_maxParallel(maxParallel)
{

}

AnalysisGraph::~AnalysisGraph()
{

}
