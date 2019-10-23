#include "testhelpers.h"

#include "analysisrunner.h"
#include "data.hpp"
#include "options.hpp"

#include <filesystem>

namespace fs = std::filesystem;

testing::AssertionResult testing::IsValidPreprocessOutPut(const Options &options)
{
    fs::path ppFile(ppFileForType(options));
    fs::path ppIndexFile(ppIndexFileForType(options));

    // Validate the output
    if (!fs::exists(ppFile))
        return testing::AssertionFailure() << "Preprocessed file does not exist:" << ppFile << std::endl;
    std::error_code ec;
    if (fs::file_size(ppFile, ec) <= 0 || ec)
        return testing::AssertionFailure() << "Preprocessed file should be greater than 0. " << ec.message() << std::endl;

    // Validate the index file
    Data data;
    AnalysisRunner::readMetaData(data, options);

    if (!fs::exists(ppIndexFile))
        return testing::AssertionFailure() << "Preprocessed index file does not exist:" << ppIndexFile << std::endl;
    const auto expectedSize = sizeof (IndexEntry) * data.numSnps;
    const auto actualSize = fs::file_size(ppIndexFile, ec);
    if (ec)
        return testing::AssertionFailure() << "Failed to get file size of "<< ppIndexFile << ": " << ec.message() << std::endl;

    if (expectedSize != actualSize)
        return testing::AssertionFailure() << "Preprocessed index file size was incorrect" << std::endl
                                           << "Expected: " << expectedSize << std::endl
                                           << "Actual: " << actualSize << std::endl;

    data.mapPreprocessBedFile(options);

    if (data.ppbedIndex.size() != data.numSnps)
        return testing::AssertionFailure() << "Preprocessed index file does not contain the correct number of entries" << std::endl
                                           << "Expected: " << data.numSnps << std::endl
                                           << "Actual: " << data.ppbedIndex.size() << std::endl;

    const auto subset = data.getMarkerIndexList();
    if (subset.empty())
        return testing::AssertionFailure() << "Marker index list should not be empty" << std::endl;

    // Test that the index entries in the subset are valid
    for (const auto marker : subset) {
        const auto index = data.ppbedIndex[marker];
        if (marker == subset.front()) {
            if (index.pos != 0)
                return testing::AssertionFailure() << "First index position should be 0" << std::endl;
        } else {
            if (index.pos <= 0)
                return testing::AssertionFailure() << "Index position should be greater than 0" << std::endl;
        }

        if (index.originalSize <= 0)
            return testing::AssertionFailure() << "Index size should be greater than 0" << std::endl;
    }

    if (subset.size() == data.numSnps)
        return testing::AssertionSuccess();

    // Test that the index entries outside the subset are valid
    auto isInvalidIndexEntry = [](const IndexEntry& index) {
        return index.pos == 0 && index.originalSize == 0;
    };

    for (unsigned int i = 0; i < subset.front(); ++i) {
        if (!isInvalidIndexEntry(data.ppbedIndex[i]))
            return testing::AssertionFailure() << "Valid index entry before the subset (should be invalid):" << std::endl
                                               << "Pos: " << data.ppbedIndex[i].pos << std::endl
                                               << "Original size: " <<  data.ppbedIndex[i].originalSize << std::endl;
    }

    for (unsigned int i = subset.back() + 1; i < data.numSnps; ++i) {
        if (!isInvalidIndexEntry(data.ppbedIndex[i]))
            return testing::AssertionFailure() << "Valid index entry after the subset (should be invalid):" << std::endl
                                               << "Pos: " << data.ppbedIndex[i].pos << std::endl
                                               << "Original size: " <<  data.ppbedIndex[i].originalSize << std::endl;
    }

    return testing::AssertionSuccess();
}
