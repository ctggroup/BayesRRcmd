#include <gtest/gtest.h>
#include <filesystem>

#include "analysisrunner.h"
#include "common.h"
#include "options.hpp"

namespace fs = std::filesystem;

class PreprocessBed : public ::testing::TestWithParam<std::tuple<PreprocessDataType, bool>> {};

TEST_P(PreprocessBed, WithAndWithoutCompressionForEachPreprocessDataType) {
    const auto params = GetParam();

    const std::string testDataDir(TEST_DATA);
    Options options;
    options.analysisType = AnalysisType::Preprocess;
    options.dataFile = testDataDir + "uk10k_chr1_1mb.bed";
    options.inputType = InputType::BED;
    options.phenotypeFile = testDataDir + "test.phen";

    // Set the test specific values
    options.preprocessDataType = std::get<0>(params);
    options.compress = std::get<1>(params);

    // Clean up old files
    fs::path ppFile(ppFileForType(options.preprocessDataType, options.dataFile));
    std::error_code ec;
    if (fs::exists(ppFile))
        ASSERT_TRUE(fs::remove(ppFile, ec)) << ec.message();

    fs::path ppIndexFile(ppIndexFileForType(options.preprocessDataType, options.dataFile));
    if (fs::exists(ppIndexFile))
        ASSERT_TRUE(fs::remove(ppIndexFile, ec)) << ec.message();

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    ASSERT_TRUE(fs::exists(ppFile));
    ASSERT_GT(fs::file_size(ppFile, ec), 0) << ec.message();

    ASSERT_TRUE(fs::exists(ppIndexFile));
    ASSERT_GT(fs::file_size(ppIndexFile, ec), 0) << ec.message();
}

INSTANTIATE_TEST_SUITE_P(PreprocessTests,
                         PreprocessBed,
                         ::testing::Combine(
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseEigen,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::Bool()));

class DISABLED_PreprocessCsvDense : public ::testing::TestWithParam<bool> {};

TEST_P(DISABLED_PreprocessCsvDense, WithAndWithoutCompression) {
    const std::string testDataDir(TEST_DATA);
    Options options;
    options.analysisType = AnalysisType::Preprocess;
    options.dataFile = testDataDir + "uk10k_chr1_1mb_transpose.csv";
    options.inputType = InputType::CSV;
    options.phenotypeFile = testDataDir + "test.csvphen";

    // Set the test specific values
    options.compress = GetParam();

    // Clean up old files
    fs::path ppFile(ppFileForType(options.preprocessDataType, options.dataFile));
    std::error_code ec;
    if (fs::exists(ppFile))
        ASSERT_TRUE(fs::remove(ppFile, ec)) << ec.message();

    fs::path ppIndexFile(ppIndexFileForType(options.preprocessDataType, options.dataFile));
    if (fs::exists(ppIndexFile))
        ASSERT_TRUE(fs::remove(ppIndexFile, ec)) << ec.message();

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    ASSERT_TRUE(fs::exists(ppFile));
    ASSERT_GT(fs::file_size(ppFile, ec), 0) << ec.message();

    ASSERT_TRUE(fs::exists(ppIndexFile));
    ASSERT_GT(fs::file_size(ppIndexFile, ec), 0) << ec.message();
}

INSTANTIATE_TEST_SUITE_P(PreprocessTests,
                         DISABLED_PreprocessCsvDense,
                         ::testing::Bool());

class PreprocessCsvSparse : public ::testing::TestWithParam<PreprocessDataType> {};

TEST_P(PreprocessCsvSparse, ExpectingFailureForSparseDataTypes) {
    const std::string testDataDir(TEST_DATA);
    Options options;
    options.analysisType = AnalysisType::Preprocess;
    options.dataFile = testDataDir + "uk10k_chr1_1mb_transpose.csv";
    options.inputType = InputType::CSV;
    options.phenotypeFile = testDataDir + "test.csvphen";

    // Set the test specific values
    options.preprocessDataType = GetParam();

    // Preprocess
    ASSERT_FALSE(AnalysisRunner::run(options)) << "Sparse preprocess data types should not be supported with CSV data";
}

INSTANTIATE_TEST_SUITE_P(PreprocessTests,
                         PreprocessCsvSparse,
                         ::testing::ValuesIn({PreprocessDataType::SparseEigen,
                                              PreprocessDataType::SparseRagged}));

class PpBayesBase : public ::testing::Test {
protected:
    Options options;

    void SetUp() override {
        const std::string testResults(TEST_RESULTS);

        options.analysisType = AnalysisType::Preprocess;
        options.chainLength = 10;
        options.burnin = 0;
        options.thin = 1;
        options.mcmcSampleFile = testResults + "ppbayes_output.csv";

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

class PpBayesBed :
        public PpBayesBase,
        public ::testing::WithParamInterface<std::tuple<AnalysisType, PreprocessDataType, bool>> {
protected:
    void SetUp() {
        PpBayesBase::SetUp();

        const std::string testDataDir(TEST_DATA);
        options.dataFile = testDataDir + "uk10k_chr1_1mb.bed";
        options.inputType = InputType::BED;
        options.phenotypeFile = testDataDir + "test.phen";
    }
};

TEST_P(PpBayesBed, SmokeTests) {

    const auto params = GetParam();
    options.preprocessDataType = std::get<1>(params);
    options.compress = std::get<2>(params);

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Run analysis
    options.analysisType = std::get<0>(params);
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    validate(options.mcmcSampleFile);
}

INSTANTIATE_TEST_SUITE_P(AnalysisSmokeTests,
                         PpBayesBed,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::PpBayes,
                                                  AnalysisType::AsyncPpBayes}),
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseEigen,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::Bool()));

class PpBayesBedGroups :
        public PpBayesBase,
        public ::testing::WithParamInterface<std::tuple<AnalysisType, PreprocessDataType, bool>> {
protected:
    void SetUp() {
        PpBayesBase::SetUp();

        const std::string testDataDir(GROUPS_TEST_DATA);
        options.dataFile = testDataDir + "uk10k_chr1_1mb.bed";
        options.inputType = InputType::BED;
        options.phenotypeFile = testDataDir + "test.phen";
        options.groupFile = testDataDir + "test.group";
        options.S = Options::parseVarianceComponents(testDataDir + "test.cva");
    }
};

TEST_P(PpBayesBedGroups, SmokeTests) {

    const auto params = GetParam();
    options.preprocessDataType = std::get<1>(params);
    options.compress = std::get<2>(params);

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Run analysis
    options.analysisType = std::get<0>(params);
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    validate(options.mcmcSampleFile);
}

INSTANTIATE_TEST_SUITE_P(PpBayesBedGroups,
                         PpBayesBed,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::PpBayes,
                                                  AnalysisType::AsyncPpBayes}),
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseEigen,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::Bool()));


class DISABLED_PpBayesCsv :
        public PpBayesBase,
        public ::testing::WithParamInterface<std::tuple<AnalysisType, bool>> {
protected:
    void SetUp() {
        PpBayesBase::SetUp();

        const std::string testDataDir(TEST_DATA);
        options.dataFile = testDataDir + "uk10k_chr1_1mb_transpose.csv";
        options.inputType = InputType::CSV;
        options.phenotypeFile = testDataDir + "test.csvphen";
    }
};

TEST_P(DISABLED_PpBayesCsv, SmokeTests) {

    const auto params = GetParam();
    options.compress = std::get<1>(params);

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Run analysis
    options.analysisType = std::get<0>(params);
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    validate(options.mcmcSampleFile);
}

INSTANTIATE_TEST_SUITE_P(AnalysisSmokeTests,
                         DISABLED_PpBayesCsv,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::PpBayes,
                                                  AnalysisType::AsyncPpBayes}),
                             ::testing::Bool()));
