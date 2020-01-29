#include <gtest/gtest.h>

#include "options.hpp"

namespace fs = std::filesystem;

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
        // --analysis-type split
        const char *argv[] = {"test", "--analysis-type", "split"};

        options.inputOptions(3, argv);
        ASSERT_EQ(AnalysisType::Unknown, options.analysisType);
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

TEST(OptionsTest, WorkingDirectory) {
    Options options;

    constexpr const char* rootTestDir = "working-directory-test";
    constexpr const char* testFile = "working-directory-test/file";
    constexpr const char* testDir = "working-directory-test/dir";
    constexpr const char* permissionsDir = "working-directory-test/permissionsDir";

    if (fs::exists(rootTestDir))
        ASSERT_TRUE(fs::remove_all(rootTestDir));

    fs::create_directory(rootTestDir);
    std::ofstream file(testFile); // create regular file
    fs::create_directory(testDir);
    fs::create_directory(permissionsDir);
    fs::permissions(permissionsDir, fs::perms::owner_read);

    {
        options.workingDirectory = fs::directory_entry();
        ASSERT_FALSE(options.validWorkingDirectory());
    }

    {
        options.workingDirectory = fs::directory_entry(testFile);
        ASSERT_FALSE(options.validWorkingDirectory());
    }

    {
        options.workingDirectory = fs::directory_entry(testDir);
        ASSERT_TRUE(options.validWorkingDirectory());
    }

    {
        options.workingDirectory = fs::directory_entry("working-directory-test/new-dir");
        ASSERT_TRUE(options.validWorkingDirectory());
    }

    {
        options.workingDirectory = fs::directory_entry(permissionsDir);
        ASSERT_FALSE(options.validWorkingDirectory());
    }

    {
        std::error_code ec;
        options.workingDirectory = fs::directory_entry("working-directory-test/permissionsDir/new-dir", ec);
        ASSERT_FALSE(options.validWorkingDirectory());
    }

    {
        options.dataFile = testFile;
        options.populateWorkingDirectory();
        ASSERT_TRUE(options.validWorkingDirectory());
    }
}

TEST(OptionsTest, WorkingDirectoryCommandLine) {
    Options options;

    constexpr const char* rootTestDir = "working-directory-test";
    constexpr const char* testFile = "working-directory-test/file";
    constexpr const char* testDir = "working-directory-test/dir";
    constexpr const char* permissionsDir = "working-directory-test/permissionsDir";

    if (fs::exists(rootTestDir))
        ASSERT_TRUE(fs::remove_all(rootTestDir));

    fs::create_directory(rootTestDir);
    std::ofstream file(testFile); // create regular file
    fs::create_directory(testDir);
    fs::create_directory(permissionsDir);
    fs::permissions(permissionsDir, fs::perms::owner_read);

    {
        const char *argv[] = {"test", "--working-directory", ""};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        ASSERT_STREQ("", options.workingDirectory.path().c_str());
    }

    {
        const char *argv[] = {"test", "--working-directory", testFile};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        ASSERT_STREQ("", options.workingDirectory.path().c_str());
    }

    {
        const char *argv[] = {"test", "--working-directory", testDir};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        ASSERT_STREQ(fs::canonical(testDir).c_str(), options.workingDirectory.path().c_str());
    }

    {
        constexpr const char* newDir = "working-directory-test/new-dir";
        const char *argv[] = {"test", "--working-directory", newDir};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        fs::path expected = fs::canonical(rootTestDir) / "new-dir";
        ASSERT_STREQ(expected.c_str(), options.workingDirectory.path().c_str());
    }

    {
        const char *argv[] = {"test", "--working-directory", permissionsDir};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        ASSERT_STREQ("", options.workingDirectory.path().c_str());
    }

    {
        const char *argv[] = {"test", "--working-directory", "working-directory-test/permissionsDir/new-dir"};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        ASSERT_STREQ("", options.workingDirectory.path().c_str());
    }

    {
        const char *argv[] = {"test", "--data-file", testFile};

        options.workingDirectory = fs::directory_entry();
        options.inputOptions(3, argv);
        ASSERT_STREQ(fs::canonical(rootTestDir).c_str(), options.workingDirectory.path().c_str());
    }
}

TEST(OptionsTest, FixedEffects) {
    Options options;
    ASSERT_TRUE(options.fixedFile.empty());
    ASSERT_EQ(0, options.fixedEffectNumber);

    {
        const char *argv[] = {"test", "--fixed_effects", "fixedEffectsFile1.csv"};

        options.inputOptions(3, argv);
        ASSERT_EQ("fixedEffectsFile1.csv", options.fixedFile);
    }

    {
        const char *argv[] = {"test", "--fixedEffectNumber", "2"};

        options.inputOptions(3, argv);
        ASSERT_EQ(2, options.fixedEffectNumber);
    }

    {
        const char *argv[] = {"test", "--fixed_effects", "fixedEffectsFile2.csv", "--fixedEffectNumber", "3"};

        options.inputOptions(5, argv);
        ASSERT_EQ("fixedEffectsFile2.csv", options.fixedFile);
        ASSERT_EQ(3, options.fixedEffectNumber);
    }
}
