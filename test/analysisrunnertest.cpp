#include <gtest/gtest.h>

#include "analysisrunner.h"
#include "limitsequencegraph.hpp"
#include "options.hpp"
#include "parallelgraph.h"
#include "sequential.h"

#include <memory>

class MakeAnalysisGraphTest : public testing::TestWithParam<std::tuple<AnalysisType, bool>> {};

TEST_P(MakeAnalysisGraphTest, makeAnalysisGraph) {
    const auto params = GetParam();
    Options options;
    options.analysisType = std::get<0>(params);
    options.useMarkerCache = std::get<1>(params);

    auto graph = AnalysisRunner::makeAnalysisGraph(options);

    switch (options.analysisType) {
    case AnalysisType::Unknown:
    {
        ASSERT_EQ(nullptr, graph);
        break;
    }
    case AnalysisType::Preprocess:
    {
        ASSERT_EQ(nullptr, graph);
        break;
    }
    case AnalysisType::PpBayes:
        // Fall through
    case AnalysisType::Gauss:
    {
        if (options.useMarkerCache) {
            auto sequential = dynamic_cast<::Sequential*>(graph.get());
            ASSERT_TRUE(sequential);
        } else {
            auto limitSequenceGraph = dynamic_cast<LimitSequenceGraph*>(graph.get());
            ASSERT_TRUE(limitSequenceGraph);
        }
        break;
    }
    case AnalysisType::AsyncPpBayes:
        // Fall through
    case AnalysisType::AsyncGauss:
    {
        auto parallelGraph = dynamic_cast<ParallelGraph*>(graph.get());
        ASSERT_TRUE(parallelGraph);
        break;
    }
    default:
        // The test does not support this analysis type
        FAIL();
    }
}

INSTANTIATE_TEST_SUITE_P(AnalysisRunnerTests, MakeAnalysisGraphTest,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::Unknown,
                                                  AnalysisType::Preprocess,
                                                  AnalysisType::PpBayes,
                                                  AnalysisType::AsyncPpBayes,
                                                  AnalysisType::Gauss,
                                                  AnalysisType::AsyncGauss}),
                             ::testing::Bool()));
