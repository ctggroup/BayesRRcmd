#include <gtest/gtest.h>

#include "analysisrunner.h"
#include "limitsequencegraph.hpp"
#include "options.hpp"
#include "parallelgraph.h"

#include <memory>

TEST(AnalysisRunnerTest, makeAnalysisGraph) {

    Options options;

    {
        // AnalysisType::Unknown
        auto graph = AnalysisRunner::makeAnalysisGraph(options);
        ASSERT_EQ(nullptr, graph);
    }

    {
        options.analysisType = AnalysisType::Preprocess;
        auto graph = AnalysisRunner::makeAnalysisGraph(options);
        ASSERT_EQ(nullptr, graph);
    }

    {
        options.analysisType = AnalysisType::PpBayes;
        auto graph = AnalysisRunner::makeAnalysisGraph(options);
        auto limitSequenceGraph = dynamic_cast<LimitSequenceGraph*>(graph.get());
        ASSERT_TRUE(limitSequenceGraph);
    }

    {
        options.analysisType = AnalysisType::AsyncPpBayes;
        auto graph = AnalysisRunner::makeAnalysisGraph(options);
        auto parallelGraph = dynamic_cast<ParallelGraph*>(graph.get());
        ASSERT_TRUE(parallelGraph);
    }

    {
        options.analysisType = AnalysisType::Gauss;
        auto graph = AnalysisRunner::makeAnalysisGraph(options);
        auto limitSequenceGraph = dynamic_cast<LimitSequenceGraph*>(graph.get());
        ASSERT_TRUE(limitSequenceGraph);
    }

    {
        options.analysisType = AnalysisType::AsyncGauss;
        auto graph = AnalysisRunner::makeAnalysisGraph(options);
        auto parallelGraph = dynamic_cast<ParallelGraph*>(graph.get());
        ASSERT_TRUE(parallelGraph);
    }
}
