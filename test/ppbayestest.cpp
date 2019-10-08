#include <gtest/gtest.h>
#include <filesystem>

#include "analysisrunner.h"
#include "common.h"
#include "options.hpp"

namespace fs = std::filesystem;

class PreprocessBed : public ::testing::TestWithParam<std::tuple<PreprocessDataType, bool, fs::directory_entry>> {};

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
    options.workingDirectory = std::get<2>(params);
    options.populateWorkingDirectory();
    ASSERT_TRUE(options.validWorkingDirectory());

    // Clean up old files
    fs::path ppFile(ppFileForType(options));
    std::error_code ec;
    if (fs::exists(ppFile))
        ASSERT_TRUE(fs::remove(ppFile, ec)) << ec.message();

    fs::path ppIndexFile(ppIndexFileForType(options));
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
                             ::testing::Bool(), // compress
                             ::testing::ValuesIn({fs::directory_entry(),
                                                  fs::directory_entry(WORKING_DIRECTORY)})));

class PreprocessCsvDense : public ::testing::TestWithParam<std::tuple<bool, fs::directory_entry>> {};

TEST_P(PreprocessCsvDense, WithAndWithoutCompression) {
    const std::string testDataDir(TEST_DATA);
    Options options;
    options.analysisType = AnalysisType::Preprocess;
    options.dataFile = testDataDir + "small_test.csv";
    options.inputType = InputType::CSV;
    options.phenotypeFile = testDataDir + "small_test.phencsv";

    // Set the test specific values
    const auto params = GetParam();
    options.compress = std::get<0>(params);
    options.workingDirectory = std::get<1>(params);
    options.populateWorkingDirectory();
    ASSERT_TRUE(options.validWorkingDirectory());

    // Clean up old files
    fs::path ppFile(ppFileForType(options));
    std::error_code ec;
    if (fs::exists(ppFile))
        ASSERT_TRUE(fs::remove(ppFile, ec)) << ec.message();

    fs::path ppIndexFile(ppIndexFileForType(options));
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
                         PreprocessCsvDense,
                         ::testing::Combine(
                             ::testing::Bool(), // compress
                             ::testing::ValuesIn({fs::directory_entry(),
                                                  fs::directory_entry(WORKING_DIRECTORY)})));

class PreprocessCsvSparse : public ::testing::TestWithParam<std::tuple<PreprocessDataType, fs::directory_entry>> {};

TEST_P(PreprocessCsvSparse, ExpectingFailureForSparseDataTypes) {
    const std::string testDataDir(TEST_DATA);
    Options options;
    options.analysisType = AnalysisType::Preprocess;
    options.dataFile = testDataDir + "uk10k_chr1_1mb_transpose.csv";
    options.populateWorkingDirectory();
    options.inputType = InputType::CSV;
    options.phenotypeFile = testDataDir + "test.csvphen";

    // Set the test specific values
    const auto params = GetParam();
    options.preprocessDataType = std::get<0>(params);
    options.workingDirectory = std::get<1>(params);
    options.populateWorkingDirectory();
    ASSERT_TRUE(options.validWorkingDirectory());

    // Preprocess
    ASSERT_FALSE(AnalysisRunner::run(options)) << "Sparse preprocess data types should not be supported with CSV data";
}

INSTANTIATE_TEST_SUITE_P(PreprocessTests,
                         PreprocessCsvSparse,
                         ::testing::Combine(
                             ::testing::ValuesIn({PreprocessDataType::SparseEigen,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::ValuesIn({fs::directory_entry(),
                                                  fs::directory_entry(WORKING_DIRECTORY)})));

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
        public ::testing::WithParamInterface<std::tuple<AnalysisType, PreprocessDataType, bool, bool, MarkerSubset>> {
protected:
    void SetUp() {
        PpBayesBase::SetUp();

        const std::string testDataDir(TEST_DATA);
        options.dataFile = testDataDir + "uk10k_chr1_1mb.bed";
        options.populateWorkingDirectory();
        options.inputType = InputType::BED;
        options.phenotypeFile = testDataDir + "test.phen";
    }
};

TEST_P(PpBayesBed, SmokeTests) {

    const auto params = GetParam();
    options.preprocessDataType = std::get<1>(params);
    options.compress = std::get<2>(params);
    options.useMarkerCache = std::get<3>(params);

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Run analysis
    options.analysisType = std::get<0>(params);
    const auto subset = std::get<4>(params);
    options.markerSubset = subset;
    std::cout << "Using marker range " << subset.first()  << " to " << subset.last() << endl;
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
                             ::testing::Bool(), // compress
                             ::testing::Bool(), // useMarkerCache
                             ::testing::ValuesIn({ MarkerSubset {0, 0},
                                                   MarkerSubset {100, 200}})
                             ));

class PpBayesBedGroups :
        public PpBayesBase,
        public ::testing::WithParamInterface<std::tuple<AnalysisType, PreprocessDataType, bool, bool>> {
protected:
    void SetUp() override {
        PpBayesBase::SetUp();

        const std::string testDataDir(GROUPS_TEST_DATA);
        options.dataFile = testDataDir + "uk10k_chr1_1mb.bed";
        options.populateWorkingDirectory();
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
    options.useMarkerCache = std::get<3>(params);

    // Preprocess
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Run analysis
    options.analysisType = std::get<0>(params);
    ASSERT_TRUE(AnalysisRunner::run(options));

    // Validate the output
    validate(options.mcmcSampleFile);
}

INSTANTIATE_TEST_SUITE_P(PpBayesBedGroups,
                         PpBayesBedGroups,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::PpBayes,
                                                  AnalysisType::AsyncPpBayes}),
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseEigen,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::Bool(), // compress
                             ::testing::Bool())); // useMarkerCache


class PpBayesCsv :
        public PpBayesBase,
        public ::testing::WithParamInterface<std::tuple<AnalysisType, bool>> {
protected:
    void SetUp() override {
        PpBayesBase::SetUp();

        const std::string testDataDir(TEST_DATA);
        options.dataFile = testDataDir + "small_test.csv";
        options.populateWorkingDirectory();
        options.inputType = InputType::CSV;
        options.phenotypeFile = testDataDir + "small_test.phencsv";
    }
};

TEST_P(PpBayesCsv, SmokeTests) {

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
                         PpBayesCsv,
                         ::testing::Combine(
                             ::testing::ValuesIn({AnalysisType::PpBayes,
                                                  AnalysisType::AsyncPpBayes}),
                             ::testing::Bool()));
