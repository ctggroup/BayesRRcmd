#include <gtest/gtest.h>

#include "options.hpp"

TEST(PreprocessedFileSplitterTest, CommandlineOptions)
{
    Options options;
    {
        // --split
        const char *argv[] = {"test", "--split"};

        options.inputOptions(2, argv);
        ASSERT_EQ(AnalysisType::Split, options.analysisType);
    }
}
