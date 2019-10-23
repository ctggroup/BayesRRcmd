#ifndef PREPROCESSEDFILESPLITTER_H
#define PREPROCESSEDFILESPLITTER_H

#include "common.h"

class Data;
class Options;
class MarkerSubset;

class PreprocessedFileSplitter
{
public:
    bool canSplit(const Options& options, const Data* data) const;
    bool split(const Options& options, const Data* data, const MarkerSubset& subset) const;

private:
    bool splitPpFile(const Options& options, const Data* data, const MarkerIndexList &markers) const;
    bool splitIndexFile(const Options& options, const Data* data, const MarkerIndexList &markers) const;
    bool writeSubsetFile(const Options& options, const Data* data, const MarkerSubset &subset) const;
};

#endif // PREPROCESSEDFILESPLITTER_H
