#include <gtest/gtest.h>
#include <filesystem>

#include "analysisrunner.h"
#include "common.h"
#include "options.hpp"
#include "testhelpers.h"

namespace fs = std::filesystem;

class BayesWBaseTest : public ::testing::Test {
protected:
    Options options;

    void SetUp() override {
        const std::string testResults(TEST_RESULTS);

        options.analysisType = AnalysisType::Preprocess;
        options.chainLength = 5;
        options.burnin = 0;
        options.thin = 1;
        options.mcmcSampleFile = testResults + "bayesw_output.csv";

        // Clean up old files
        fs::path resultsFile(options.mcmcSampleFile);
        std::error_code ec;
        if (fs::exists(resultsFile))
            ASSERT_TRUE(fs::remove(resultsFile, ec)) << ec.message();

        if (!fs::exists(resultsFile.parent_path()))
            ASSERT_TRUE(fs::create_directories(resultsFile.parent_path(), ec)) << ec.message();
    }

    virtual void validate(const fs::path &resultsFile) {
        ASSERT_TRUE(fs::exists(resultsFile));
        std::error_code ec;
        ASSERT_GT(fs::file_size(resultsFile, ec), 0) << ec.message();

        std::ifstream stream(resultsFile);
        const auto lineCount = std::count(std::istreambuf_iterator<char>(stream),
                                          std::istreambuf_iterator<char>(), '\n');
        const auto expectedLineCount = 1 + std::ceil(static_cast<double>(options.chainLength) /
                                                     static_cast<double>(options.thin));
        ASSERT_EQ(lineCount, expectedLineCount);
    }
};

using BayesWTestParams = std::tuple<AnalysisType, PreprocessDataType, bool, bool, testing::FixedEffectsParams>;

class BayesWTest :
        public BayesWBaseTest,
        public ::testing::WithParamInterface<BayesWTestParams> {
protected:
    void SetUp() {
        BayesWBaseTest::SetUp();

        const std::string testDataDir(GAUSS_TEST_DATA);
        options.dataFile = testDataDir + "data.bed";
        options.populateWorkingDirectory();
        options.inputType = InputType::BED;
        options.failureFile = testDataDir + "data.fail";
        options.quad_points = "7";
        options.S = MatrixXd(1, 2);
        options.S << 0.01, 0.1;
        options.phenotypeFile = testDataDir + "data.phen";
    }
};

TEST_P(BayesWTest, SmokeTests) {

    const auto params = GetParam();
    options.preprocessDataType = std::get<1>(params);
    options.compress = std::get<2>(params);
    options.useMarkerCache = std::get<3>(params);

    const auto fixedEffectsParams = std::get<4>(params);
    if (fixedEffectsParams.fixedEffectNumber > 0) {
        options.fixedEffectNumber = fixedEffectsParams.fixedEffectNumber;
        const std::string testDataDir(GAUSS_TEST_DATA);
        options.fixedFile = testDataDir + fixedEffectsParams.fixedFile;
    }

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Run analysis
    options.analysisType = std::get<0>(params);
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    validate(options.mcmcSampleFile);
}

INSTANTIATE_TEST_SUITE_P(AnalysisSmokeTests,
                         BayesWTest,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::Gauss,
                                                  AnalysisType::AsyncGauss}),
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::Bool(), // compress
                             ::testing::Bool(), // useMarkerCache
                             ::testing::ValuesIn({testing::FixedEffectsParams{0, ""},
                                                  testing::FixedEffectsParams{1, "1FE-1000IND-1NA.csv"}}
                                                 )));
