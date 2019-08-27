#ifndef ANALYSISRUNNER_H
#define ANALYSISRUNNER_H

#include <memory>

struct Options;

class AnalysisGraph;

namespace AnalysisRunner {
std::unique_ptr<AnalysisGraph> makeAnalysisGraph(const Options &options);
bool run(const Options &options);
}

#endif // ANALYSISRUNNER_H
