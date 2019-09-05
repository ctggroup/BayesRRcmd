#include <gtest/gtest.h>

#include "options.hpp"
#include "data.hpp"

TEST(OptionsTest, VarianceFromCommandLine) {
    Options options;

    {
        // Default values
        MatrixXd expected(1, 3);
        expected << 0.01, 0.001, 0.0001;

        ASSERT_EQ(expected, options.S);
    }

    {
        // 1 group, 2 components
        MatrixXd expected(1, 2);
        expected << 0.1, 0.1;

        const char *argv[] = {"test", "--S", "0.1,0.1"};

        options.inputOptions(3, argv);
        ASSERT_EQ(expected, options.S);
    }

    {
        // 2 groups, 3 components
        MatrixXd expected(2, 3);
        expected << 0.01, 0.001, 0.0001, 0.02, 0.002, 0.0002;

        const char *argv[] = {"test", "--S", "0.01,0.001,0.0001;0.02,0.002,0.0002"};

        options.inputOptions(3, argv);
        ASSERT_EQ(expected, options.S);
    }

    {
        // From file - 2 groups, 2 components
        MatrixXd expected(2, 2);
        expected << 0.001, 0.01, 0.001, 0.01;

        std::string testFile {GROUPS_TEST_DATA};
        testFile += "test.cva";
        {
            ifstream check(testFile);
            ASSERT_TRUE(check.is_open()) << "Failed to open: " << testFile;
        }

        const char *argv[] = {"test", "--S", testFile.c_str()};

        options.inputOptions(3, argv);
        ASSERT_EQ(expected, options.S);
    }

    {
        // File does not exist
        // Default constructed MatrixXd is the expected value for an error
        MatrixXd expected;

        std::string testFile {GROUPS_TEST_DATA};
        testFile += "non-existant.cva";
        {
            ifstream check(testFile);
            ASSERT_FALSE(check.is_open()) << "File should not exist: " << testFile;
        }

        const char *argv[] = {"test", "--S", testFile.c_str()};

        options.inputOptions(3, argv);
        ASSERT_EQ(expected, options.S);
    }
}

TEST(OptionsTest, DataFileType) {
    Options options;
    ASSERT_EQ(InputType::Unknown, options.inputType);

    {
        // BED
        const char *argv[] = {"test", "--data-file", "foo.bed"};

        options.inputOptions(3, argv);
        ASSERT_EQ(InputType::BED, options.inputType);
    }

    {
        // CSV
        const char *argv[] = {"test", "--data-file", "foo.csv"};

        options.inputOptions(3, argv);
        ASSERT_EQ(InputType::CSV, options.inputType);
    }

    {
        // BED
        const char *argv[] = {"test", "--data-file", "foo"};

        options.inputOptions(3, argv);
        ASSERT_EQ(InputType::Unknown, options.inputType);
    }
}

TEST(OptionsTest, AnalysisType) {
    Options options;
    ASSERT_EQ(AnalysisType::Unknown, options.analysisType);

    {
        // --preprocess
        const char *argv[] = {"test", "--preprocess"};

        options.inputOptions(2, argv);
        ASSERT_EQ(AnalysisType::Preprocess, options.analysisType);
    }

    {
        // --analysis-type preprocess
        const char *argv[] = {"test", "--analysis-type", "preprocess"};

        options.inputOptions(3, argv);
        ASSERT_EQ(AnalysisType::Preprocess, options.analysisType);
    }

    {
        // --analysis-type ppbayes
        const char *argv[] = {"test", "--analysis-type", "ppbayes"};

        options.inputOptions(3, argv);
        ASSERT_EQ(AnalysisType::PpBayes, options.analysisType);
    }

    {
        // --analysis-type asyncppbayes
        const char *argv[] = {"test", "--analysis-type", "asyncppbayes"};

        options.inputOptions(3, argv);
        ASSERT_EQ(AnalysisType::AsyncPpBayes, options.analysisType);
    }

    {
        // --analysis-type gauss
        const char *argv[] = {"test", "--analysis-type", "gauss"};

        options.inputOptions(3, argv);
        ASSERT_EQ(AnalysisType::Gauss, options.analysisType);
    }

    {
        // --analysis-type gauss
        const char *argv[] = {"test", "--analysis-type", "asyncgauss"};

        options.inputOptions(3, argv);
        ASSERT_EQ(AnalysisType::AsyncGauss, options.analysisType);
    }
}

TEST(OptionsTest, PreprocessDataType) {
    Options options;
    ASSERT_EQ(PreprocessDataType::Dense, options.preprocessDataType);

    {
        // SparseEigen
        const char *argv[] = {"test", "--sparse-data", "eigen"};

        options.inputOptions(3, argv);
        ASSERT_EQ(PreprocessDataType::SparseEigen, options.preprocessDataType);
    }

    {
        // SparseRagged
        const char *argv[] = {"test", "--sparse-data", "ragged"};

        options.inputOptions(3, argv);
        ASSERT_EQ(PreprocessDataType::SparseRagged, options.preprocessDataType);
    }

    {
        // None
        const char *argv[] = {"test", "--sparse-data", "foo"};

        options.inputOptions(3, argv);
        ASSERT_EQ(PreprocessDataType::None, options.preprocessDataType);
    }
}

TEST(OptionsTest, MarkerSubset) {
    Data data;
    data.numSnps = 100;

    Options options;
    {
        // default markerSubset
        ASSERT_TRUE(options.validMarkerSubset(&data));

        const auto markers = options.getMarkerSubset(&data);
        ASSERT_EQ(0, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        options.markerSubset = {0, 100};
        ASSERT_TRUE(options.validMarkerSubset(&data));

        const auto markers = options.getMarkerSubset(&data);
        ASSERT_EQ(0, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        options.markerSubset = {99, 1};
        ASSERT_TRUE(options.validMarkerSubset(&data));

        const auto markers = options.getMarkerSubset(&data);
        ASSERT_EQ(99, markers.front());
        ASSERT_EQ(99, markers.back());
    }

    {
        options.markerSubset = {0, 101};
        ASSERT_FALSE(options.validMarkerSubset(&data));
    }

    {
        options.markerSubset = {99, 2};
        ASSERT_FALSE(options.validMarkerSubset(&data));
    }

    {
        options.markerSubset = {100, 1};
        ASSERT_FALSE(options.validMarkerSubset(&data));
    }

    {
        options.markerSubset = {25, 50};
        ASSERT_TRUE(options.validMarkerSubset(&data));

        const auto markers = options.getMarkerSubset(&data);
        ASSERT_EQ(25, markers.front());
        ASSERT_EQ(74, markers.back());
    }
}
