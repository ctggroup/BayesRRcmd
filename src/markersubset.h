#ifndef MARKERSUBSET_H
#define MARKERSUBSET_H

#include "common.h"

struct MarkerSubset {
    MarkerSubset(unsigned int start = 0, unsigned int size = 0)
        : m_start(start)
        , m_size(size)
    {}

    auto asTuple() const { return std::tie(m_start, m_size); }

    unsigned int first() const { return m_start; }
    unsigned long last() const;
    unsigned int size() const { return m_size; }

    void clamp(unsigned int markerCount);

    bool isValid(unsigned int markerCount) const;
    MarkerIndexList toMarkerIndexList(unsigned int markerCount) const;

private:
    unsigned int m_start = 0;
    unsigned int m_size = 0;
};

bool operator==(const MarkerSubset &lhs, const MarkerSubset &rhs);
bool operator!=(const MarkerSubset &lhs, const MarkerSubset &rhs);

std::ostream& operator<<(std::ostream& out, const MarkerSubset& subset);
std::istream& operator>>(std::istream &in,  MarkerSubset& subset);

// Checks that each marker is sampled exactly once
bool isValid(const std::vector<MarkerSubset> subsets, unsigned int markerCount);

std::vector<MarkerSubset> generateEqualSubsets(unsigned int subsetCount,
                                               unsigned int markerCount);

#endif // MARKERSUBSET_H
