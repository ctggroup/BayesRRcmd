#include <gtest/gtest.h>

#include "data.hpp"

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

    data.numInds = individualCount;
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
                             std::make_tuple(10, 5, 3, "fixed-effects/5FE-10IND-3NA.csv")
                             ));
