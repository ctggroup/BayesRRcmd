#ifndef MARKERSUBSET_H
#define MARKERSUBSET_H

#include "common.h"

struct MarkerSubset {
    unsigned int start = 0;
    unsigned int size = 0;

    auto asTuple() const { return std::tie(start, size); }

    unsigned int first() const { return start; }
    unsigned int last() const { return start + size - 1; }

    void clamp(unsigned int markerCount);

    bool isValid(unsigned int markerCount) const;
    MarkerIndexList toMarkerIndexList(unsigned int markerCount) const;
};

bool operator==(const MarkerSubset &lhs, const MarkerSubset &rhs);
bool operator!=(const MarkerSubset &lhs, const MarkerSubset &rhs);

// Checks that each marker is sampled exactly once
bool isValid(const std::vector<MarkerSubset> subsets, unsigned int markerCount);

std::vector<MarkerSubset> generateEqualSubsets(unsigned int subsetCount,
                                               unsigned int markerCount);

#endif // MARKERSUBSET_H
