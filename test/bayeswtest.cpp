#include <gtest/gtest.h>
#include <filesystem>

#include "analysisrunner.h"
#include "common.h"
#include "options.hpp"

namespace fs = std::filesystem;

class BayesWBaseTest : public ::testing::Test {
protected:
    Options options;

    void SetUp() override {
        const std::string testResults(TEST_RESULTS);

        options.chainLength = 4; // BayesW is super slow - increase after performance improvements
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

class BayesWTest :
        public BayesWBaseTest,
        public ::testing::WithParamInterface<std::tuple<PreprocessDataType>> {
protected:
    void SetUp() {
        BayesWBaseTest::SetUp();

        const std::string testDataDir(GAUSS_TEST_DATA);
        options.analysisType = AnalysisType::Gauss;
        options.dataFile = testDataDir + "data.bed";
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
    options.preprocessDataType = std::get<0>(params);

    // Run analysis
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    validate(options.mcmcSampleFile);
}

INSTANTIATE_TEST_SUITE_P(AnalysisSmokeTests,
                         BayesWTest,
                         ::testing::Combine(
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseRagged})));
