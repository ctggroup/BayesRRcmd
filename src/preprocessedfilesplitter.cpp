#include "preprocessedfilesplitter.h"

#include "data.hpp"
#include "options.hpp"

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <iostream>

namespace fs = std::filesystem;

bool PreprocessedFileSplitter::canSplit(const Options &options, const Data *data) const
{
    fs::path ppFilePath = ppFileForType(options);
    if (!fs::exists(ppFilePath)) {
        std::cerr << "Source file does not exist: " << ppFilePath << std::endl;
        return false;
    }

    fs::path outPath = ppFileForType(options, options.splitDestination);
    if (!validateDirectory(outPath.parent_path())) {
        std::cerr << "Invalid split destination:" << outPath << std::endl;
        return false;
    }

    std::error_code ec;
    ppFilePath = fs::canonical(ppFilePath.parent_path(), ec);
    if (ec) {
        std::cerr << ec.message() << std::endl;
        return false;
    }

    outPath = fs::canonical(outPath.parent_path(), ec);
    if (ec) {
        std::cerr << ec.message() << std::endl;
        return false;
    }

    if (outPath == ppFilePath) {
        std::cerr << "Destination directory is the same as the input data directory path. " << std::endl
                  << "Source" << ppFilePath << std::endl
                  << "Destination: " << outPath << std::endl;
        return false;
    }

    return true;
}

bool PreprocessedFileSplitter::split(const Options &options, const Data *data, const MarkerSubset &subset) const
{
    if (!canSplit(options, data))
        return false;

    // auto markers = subset.toMarkerIndexList(data->numSnps);
    auto markers = data->getMarkerIndexList();
    if (markers.empty()) {
        std::cerr << "Invalid marker subset:" << subset << std::endl;
        return false;
    }

    int rank = 0;
#ifdef MPI_ENABLED
    if (options.useHybridMpi)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    std::cout << rank << " Splitting preprocessed data:" << std::endl
              << " - source: " << ppFileForType(options) << std::endl
              << " - destination: " << ppFileForType(options, options.splitDestination) << std::endl
              << " - markers: " << markers.front() << " to " << markers.back() << " of " << data->numSnps << std::endl;

    if (!splitPpFile(options, data, markers))
        return false;

    if (!splitIndexFile(options, data, markers))
        return false;

    if (!writeSubsetFile(options, data, subset))
        return false;

    return true;
}

bool PreprocessedFileSplitter::splitPpFile(const Options &options, const Data *data, const MarkerIndexList &markers) const
{
    const fs::path source = ppFileForType(options);
    const fs::path outPath = ppFileForType(options, options.splitDestination);

    auto output = std::ofstream(outPath.c_str(), ios::binary);
    if (output.fail()){
        cerr << "Error: Unable to open the output file for writing: " << outPath  << endl;
        return false;
    }

    for (const auto i : markers) {
        const auto index = data->ppbedIndex[i];
        auto *input = reinterpret_cast<char*>(data->ppBedMap);
        input += index.pos;

        output.write(input, index.compressedSize);

        if (output.fail()) {
            std::cerr << "Error whilst splitting preprocessed data" << std::endl;
            break;
        }
    }

    output.close();
    return !output.fail();
}

bool PreprocessedFileSplitter::splitIndexFile(const Options &options, const Data *data, const MarkerIndexList &markers) const
{
    const fs::path source = ppIndexFileForType(options);
    const fs::path outPath = ppIndexFileForType(options, options.splitDestination);

    auto output = std::ofstream(outPath.c_str(), ios::binary);
    if (output.fail()){
        cerr << "Error: Unable to open the output index file for writing: " << outPath  << endl;
        return false;
    }

    static const IndexEntry invalidEntry;
    auto writeMarker = [&output] (const IndexEntry &entry) {
        output.write(reinterpret_cast<const char *>(&entry),
                             sizeof(entry));
    };

    // Fill the front of the marker index
    for (unsigned int i = 0; i < markers.front(); ++i) {
        writeMarker(invalidEntry);
    }

    // Fill the subset of the marker index
    size_t pos = 0;
    for (const auto i : markers) {
        IndexEntry entry = data->ppbedIndex[i];
        entry.pos = pos;
        writeMarker(entry);
        pos += entry.compressedSize;
    }

    // Fill the back of the marker index
    for (unsigned int i = markers.back() + 1; i < data->numSnps; ++i) {
        writeMarker(invalidEntry);
    }

    output.close();
    return !output.fail();
}

bool PreprocessedFileSplitter::writeSubsetFile(const Options &options, const Data *data, const MarkerSubset &subset) const
{
    const auto ppSubsetFile = ppSubsetFileForType(options, options.splitDestination);

    // Write the MarkerSubset this split contains
    std::ofstream output = std::ofstream(ppSubsetFile.c_str(), ios::binary);
    if (output.fail()) {
        cerr << "Error: Unable to open the preprocessed bed subset file for writing: " << ppSubsetFile << endl;
        return false;
    }

    output << subset;
    output.close();
    return !output.fail();
}
