#include <gtest/gtest.h>

#include "markersubset.h"
#include "data.hpp"

TEST(MarkerSubset, MarkerSubset) {
    MarkerSubset s = {0, 2};
    ASSERT_EQ(0, s.first());
    ASSERT_EQ(1, s.last());
    ASSERT_EQ(2, s.size());
}

TEST(MarkerSubset, isValid) {
    MarkerSubset s;
    {
        // default markerSubset
        ASSERT_TRUE(s.isValid(100));

        const auto markers = s.toMarkerIndexList(100);
        ASSERT_EQ(0, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        s = {0, 100};
        ASSERT_TRUE(s.isValid(100));

        const auto markers = s.toMarkerIndexList(100);
        ASSERT_EQ(0, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        s = {99, 1};
        ASSERT_TRUE(s.isValid(100));

        const auto markers = s.toMarkerIndexList(100);
        ASSERT_EQ(99, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        s = {25, 50};
        ASSERT_TRUE(s.isValid(100));

        const auto markers = s.toMarkerIndexList(100);
        ASSERT_EQ(25, markers.front());
        ASSERT_EQ(74, markers.back());
    }

    {
        s = {0, 101};
        ASSERT_FALSE(s.isValid(100));
    }

    {
        s = {99, 2};
        ASSERT_FALSE(s.isValid(100));
    }

    {
        s = {100, 1};
        ASSERT_FALSE(s.isValid(100));
    }
}

TEST(MarkerSubset, isValidSubsetList) {
    std::vector<MarkerSubset> subsets;
    ASSERT_FALSE(isValid(subsets, 0));

    {
        subsets = {{0, 50}, {50, 50}};
        ASSERT_TRUE(isValid(subsets, 100));
    }

    {
        // Subsets overlap
        subsets = {{0, 50}, {49, 50}};
        ASSERT_FALSE(isValid(subsets, 100));
    }

    {
        // Subsets are not contiguous
        subsets = {{0, 50}, {51, 50}};
        ASSERT_FALSE(isValid(subsets, 100));
    }

    {
        // Too few markers
        subsets = {{0, 25}, {25, 25}};
        ASSERT_FALSE(isValid(subsets, 100));
    }

    {
        // Too many markers
        subsets = {{0, 75}, {75, 75}};
        ASSERT_FALSE(isValid(subsets, 100));
    }
}

TEST(MarkerSubset, getMarkerSubset) {
    Data data;
    data.numSnps = 100;

    {
        // default markerSubset
        ASSERT_TRUE(data.validMarkerSubset());

        const auto markers = data.getMarkerIndexList();
        ASSERT_EQ(0, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        ASSERT_TRUE(data.setMarkerSubset({0, 100}));
        ASSERT_TRUE(data.validMarkerSubset());

        const auto markers = data.getMarkerIndexList();
        ASSERT_EQ(0, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        ASSERT_TRUE(data.setMarkerSubset({99, 1}));
        ASSERT_TRUE(data.validMarkerSubset());

        const auto markers = data.getMarkerIndexList();
        ASSERT_EQ(99, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        ASSERT_FALSE(data.setMarkerSubset({0, 101}));
    }

    {
        ASSERT_FALSE(data.setMarkerSubset({99, 2}));
    }

    {
        ASSERT_FALSE(data.setMarkerSubset({100, 1}));
    }

    {
        ASSERT_TRUE(data.setMarkerSubset({25, 50}));
        ASSERT_TRUE(data.validMarkerSubset());

        const auto markers = data.getMarkerIndexList();
        ASSERT_EQ(25, markers.front());
        ASSERT_EQ(74, markers.back());
    }
}

class MarkerSubsetLists : public testing::TestWithParam<std::tuple<int, int>> {};

TEST_P(MarkerSubsetLists, generateEqualSubsets) {
    const std::vector<MarkerSubset> emptyList;
    const auto params = GetParam();

    const auto markerCount = static_cast<unsigned int>(std::get<0>(params));
    const auto subsetCount = static_cast<unsigned int>(std::get<1>(params));
    const auto subsets = generateEqualSubsets(subsetCount, markerCount);

    ASSERT_NE(subsets, emptyList);
    ASSERT_EQ(subsetCount, subsets.size());
    ASSERT_TRUE(isValid(subsets, markerCount));
}

INSTANTIATE_TEST_SUITE_P(MarkerSubset,
                         MarkerSubsetLists,
                         ::testing::Combine(
                             ::testing::ValuesIn({5, 10, 97, 100, 6717}), // markerCount
                             ::testing::Range(1, 10))); // subsetCount
