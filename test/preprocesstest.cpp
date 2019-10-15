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
