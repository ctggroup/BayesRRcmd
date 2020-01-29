#include <gtest/gtest.h>

#include "data.hpp"

TEST(FixedEffectsTest, FE1NA0) {
    const std::string testDataDir = {TEST_DATA};
    const unsigned int individualCount = 10;
    const unsigned int fixedEffectCount = 1;

    Data data;

    ASSERT_EQ(0, data.X.rows());
    ASSERT_EQ(0, data.X.cols());

    data.numInds = individualCount;
    data.readCSV(testDataDir + "fixed-effects/1FE-10IND-0NA.csv", fixedEffectCount);

    ASSERT_EQ(individualCount, data.X.rows());
    ASSERT_EQ(fixedEffectCount, data.X.cols());
}

TEST(FixedEffectsTest, FE5NA0) {
    const std::string testDataDir = {TEST_DATA};
    const unsigned int individualCount = 10;
    const unsigned int fixedEffectCount = 5;

    Data data;

    ASSERT_EQ(0, data.X.rows());
    ASSERT_EQ(0, data.X.cols());

    data.numInds = individualCount;
    data.readCSV(testDataDir + "fixed-effects/5FE-10IND-0NA.csv", fixedEffectCount);

    ASSERT_EQ(individualCount, data.X.rows());
    ASSERT_EQ(fixedEffectCount, data.X.cols());
}
