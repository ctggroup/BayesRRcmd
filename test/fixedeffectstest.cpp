#include <gtest/gtest.h>

#include "analysisrunner.h"
#include "data.hpp"
#include "options.hpp"

// individual count, fixed effect count, NA count, fixed effects file
using ReadFixedEffectsParams = std::tuple<unsigned int, unsigned int, unsigned int, std::string>;

class ReadFixedEffectsParamTest : public ::testing::TestWithParam<ReadFixedEffectsParams> {
protected:
    Data data;
};

TEST_P(ReadFixedEffectsParamTest, ReadFixedEffects) {
    const auto params = GetParam();

    const unsigned int individualCount = std::get<0>(params);
    const unsigned int fixedEffectCount = std::get<1>(params);
    const unsigned int naCount = std::get<2>(params);

    const std::string testDataDir = {TEST_DATA};
    const unsigned int expectedIndividuals = individualCount - naCount;

    ASSERT_EQ(0, data.X.rows());
    ASSERT_EQ(0, data.X.cols());

    data.totalIndividuals = individualCount;
    data.readFixedEffects(testDataDir + std::get<3>(params), fixedEffectCount);

    ASSERT_EQ(expectedIndividuals, data.X.rows());
    ASSERT_EQ(fixedEffectCount, data.X.cols());

    ASSERT_EQ(individualCount * fixedEffectCount, data.m_fixedEffectsData->capacity());
    ASSERT_EQ(expectedIndividuals * fixedEffectCount, data.m_fixedEffectsData->size());
}

INSTANTIATE_TEST_SUITE_P(FixedEffectsTests,
                         ReadFixedEffectsParamTest,
                         ::testing::Values(
                             std::make_tuple(10, 1, 0, "fixed-effects/1FE-10IND-0NA.csv"),
                             std::make_tuple(10, 5, 0, "fixed-effects/5FE-10IND-0NA.csv"),
                             std::make_tuple(10, 1, 1, "fixed-effects/1FE-10IND-1NA.csv"),
                             std::make_tuple(10, 5, 3, "fixed-effects/5FE-10IND-3NA.csv"),
                             std::make_tuple(500, 5, 3, "5FE-500IND-3NA.csv"),
                             std::make_tuple(3642, 1, 0, "1FE-3642IND-0NA.csv"),
                             std::make_tuple(3642, 1, 1, "1FE-3642IND-1NA.csv"),
                             std::make_tuple(3642, 1, 3, "1FE-3642IND-3NA.csv"),
                             std::make_tuple(3642, 5, 0, "5FE-3642IND-0NA.csv"),
                             std::make_tuple(3642, 5, 1, "5FE-3642IND-1NA.csv"),
                             std::make_tuple(3642, 5, 3, "5FE-3642IND-3NA.csv")
                             ));

// InputType, data file, phenotype file, fixed effects file, fixed effects number
using ReadMetaDataParams = std::tuple<InputType, std::string, std::string, std::string, int>;

class ReadMetaDataParamTest : public ::testing::TestWithParam<ReadMetaDataParams> {
protected:
    Data data;
    Options options;
};

TEST_P(ReadMetaDataParamTest, ReadMetaData) {
    const auto params = GetParam();
    const std::string testDataDir(TEST_DATA);

    options.inputType = std::get<0>(params);

    options.dataFile = testDataDir + std::get<1>(params);
    options.phenotypeFile = testDataDir + std::get<2>(params);

    options.fixedFile = testDataDir + std::get<3>(params);
    options.fixedEffectNumber = std::get<4>(params);

    AnalysisRunner::readMetaData(data, options);

    ASSERT_EQ(data.totalIndividuals, data.indInfoVec.size());
    if (options.inputType == InputType::BED)
        ASSERT_EQ(data.totalIndividuals, data.indInfoMap.size());

    ASSERT_EQ(data.X.rows(), data.activeIndividuals);
    ASSERT_EQ(data.X.rows(), data.y.size());
}

INSTANTIATE_TEST_SUITE_P(ReadMetaDataTests,
                         ReadMetaDataParamTest,
                         ::testing::Values(
                             std::make_tuple(InputType::BED, "uk10k_chr1_1mb.bed", "test.phen", "1FE-3642IND-0NA.csv", 1),
                             std::make_tuple(InputType::BED, "uk10k_chr1_1mb.bed", "test.phen", "1FE-3642IND-1NA.csv", 1),
                             std::make_tuple(InputType::BED, "uk10k_chr1_1mb.bed", "test.phen", "1FE-3642IND-3NA.csv", 1),
                             std::make_tuple(InputType::BED, "uk10k_chr1_1mb.bed", "test.phen", "5FE-3642IND-0NA.csv", 5),
                             std::make_tuple(InputType::BED, "uk10k_chr1_1mb.bed", "test.phen", "5FE-3642IND-1NA.csv", 5),
                             std::make_tuple(InputType::BED, "uk10k_chr1_1mb.bed", "test.phen", "5FE-3642IND-3NA.csv", 5),
                             std::make_tuple(InputType::CSV, "small_test.csv", "small_test.phencsv", "5FE-500IND-3NA.csv", 5),
                             std::make_tuple(InputType::CSV, "uk10k_chr1_1mb_transpose.csv", "test.csvphen", "1FE-3642IND-0NA.csv", 1),
                             std::make_tuple(InputType::CSV, "uk10k_chr1_1mb_transpose.csv", "test.csvphen", "1FE-3642IND-1NA.csv", 1),
                             std::make_tuple(InputType::CSV, "uk10k_chr1_1mb_transpose.csv", "test.csvphen", "1FE-3642IND-3NA.csv", 1),
                             std::make_tuple(InputType::CSV, "uk10k_chr1_1mb_transpose.csv", "test.csvphen", "5FE-3642IND-0NA.csv", 5),
                             std::make_tuple(InputType::CSV, "uk10k_chr1_1mb_transpose.csv", "test.csvphen", "5FE-3642IND-1NA.csv", 5),
                             std::make_tuple(InputType::CSV, "uk10k_chr1_1mb_transpose.csv", "test.csvphen", "5FE-3642IND-3NA.csv", 5)
                             ));
