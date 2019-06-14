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
