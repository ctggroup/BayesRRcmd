#include <gtest/gtest.h>

#include "analysisrunner.h"
#include "data.hpp"
#include "options.hpp"
#include "preprocessedfilesplitter.h"
#include "testhelpers.h"

#include <filesystem>

namespace fs = std::filesystem;

TEST(PreprocessedFileSplitterTest, CommandlineOptions)
{
    Options options;
    {
        // --split
        const auto destination = "destination";
        const char *argv[] = {"test", "--split", destination};

        options.inputOptions(2, argv);
        ASSERT_EQ(AnalysisType::Split, options.analysisType);
        ASSERT_EQ(options.splitDestination, destination);
    }
}

class PreprocessedFileSplitterTest : public ::testing::TestWithParam<std::tuple<bool, fs::directory_entry, std::string, MarkerSubset>> {};

TEST_P(PreprocessedFileSplitterTest, SplitPreprocessedFile) {
    const auto params = GetParam();

    const std::string testDataDir(TEST_DATA);
    Options options;
    options.analysisType = AnalysisType::Preprocess;
    options.dataFile = testDataDir + "uk10k_chr1_1mb.bed";
    options.inputType = InputType::BED;
    options.phenotypeFile = testDataDir + "test.phen";

    // Set the test specific values
    options.preprocessDataType = PreprocessDataType::SparseRagged;
    options.compress = std::get<0>(params);
    const auto workingDirectory = std::get<1>(params);
    options.workingDirectory = workingDirectory;
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
    ASSERT_TRUE(testing::IsValidPreprocessOutPut(options));

    // Split
    options.analysisType = AnalysisType::Split;
    const std::string destination = std::get<2>(params);
    options.splitDestination = destination;
    options.preprocessSubset = std::get<3>(params);

    if (destination.empty() || (destination == workingDirectory.path())) {
        ASSERT_FALSE(AnalysisRunner::run(options));
        return;
    } else {
        ASSERT_TRUE(AnalysisRunner::run(options));
    }

    options.workingDirectory = fs::directory_entry(destination);
    options.populateWorkingDirectory();
    ASSERT_TRUE(testing::IsValidPreprocessOutPut(options));
}

INSTANTIATE_TEST_SUITE_P(PreprocessedFileSplitterTests,
                         PreprocessedFileSplitterTest,
                         ::testing::Combine(
                             ::testing::Bool(), // compress
                             ::testing::ValuesIn({fs::directory_entry(),
                                                  fs::directory_entry(WORKING_DIRECTORY)}),
                             ::testing::ValuesIn({std::string(), // invalid
                                                  std::string(WORKING_DIRECTORY), // invalid
                                                  std::string(SPLIT_DESTINATION)}), // valid
                             ::testing::ValuesIn({MarkerSubset{0, 0}, MarkerSubset{100, 100}})));
