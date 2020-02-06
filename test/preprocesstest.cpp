#include <gtest/gtest.h>
#include <filesystem>

#include "analysisrunner.h"
#include "common.h"
#include "data.hpp"
#include "options.hpp"
#include "testhelpers.h"

namespace fs = std::filesystem;

// preprocess data type, compress, working directory, marker subset, preprocess chunks, fixed effects params
using PreprocessBedParams = std::tuple<PreprocessDataType, bool, fs::directory_entry, MarkerSubset, int, testing::FixedEffectsParams>;

class PreprocessBed : public ::testing::TestWithParam<PreprocessBedParams> {};

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
    options.preprocessSubset = std::get<3>(params);
    options.preprocessChunks = std::get<4>(params);

    const auto fixedEffectsParams = std::get<5>(params);
    if (fixedEffectsParams.fixedEffectNumber > 0) {
        options.fixedEffectNumber = fixedEffectsParams.fixedEffectNumber;
        options.fixedFile = testDataDir + fixedEffectsParams.fixedFile;
    }

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
    ASSERT_TRUE(testing::IsValidPreprocessOutPut(options));
}

INSTANTIATE_TEST_SUITE_P(PreprocessTests,
                         PreprocessBed,
                         ::testing::Combine(
                             ::testing::ValuesIn({PreprocessDataType::Dense,
                                                  PreprocessDataType::SparseEigen,
                                                  PreprocessDataType::SparseRagged}),
                             ::testing::Bool(), // compress
                             ::testing::ValuesIn({fs::directory_entry(),
                                                  fs::directory_entry(WORKING_DIRECTORY)}),
                             ::testing::ValuesIn({MarkerSubset{0, 0}, MarkerSubset{100, 100}}),
                             ::testing::ValuesIn({1, 3, 7}),
                             ::testing::ValuesIn({testing::FixedEffectsParams{0, ""},
                                                  testing::FixedEffectsParams{5, "5FE-3642IND-3NA.csv"}}
                                                 )));

using PreprocessCsvDenseParams = std::tuple<bool, fs::directory_entry, testing::FixedEffectsParams>;

class PreprocessCsvDense : public ::testing::TestWithParam<PreprocessCsvDenseParams> {};

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

    const auto fixedEffectsParams = std::get<2>(params);
    if (fixedEffectsParams.fixedEffectNumber > 0) {
        options.fixedEffectNumber = fixedEffectsParams.fixedEffectNumber;
        options.fixedFile = testDataDir + fixedEffectsParams.fixedFile;
    }

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
    ASSERT_TRUE(testing::IsValidPreprocessOutPut(options));
}

INSTANTIATE_TEST_SUITE_P(PreprocessTests,
                         PreprocessCsvDense,
                         ::testing::Combine(
                             ::testing::Bool(), // compress
                             ::testing::ValuesIn({fs::directory_entry(),
                                                  fs::directory_entry(WORKING_DIRECTORY)}),
                             ::testing::ValuesIn({testing::FixedEffectsParams{0, ""},
                                                  testing::FixedEffectsParams{5, "5FE-500IND-3NA.csv"}}
                                                 )));

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
