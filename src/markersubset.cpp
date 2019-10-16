#include "markersubset.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

namespace  {
    static constexpr auto k_markerSubsetSize = sizeof (MarkerSubset);
}

unsigned long MarkerSubset::last() const
{
     if (start == 0 && size == 0)
         return 0;

     return start + size - 1;
}

void MarkerSubset::clamp(unsigned int markerCount)
{
    if (last() < markerCount)
        return;

    size = markerCount - start;
}

bool MarkerSubset::isValid(unsigned int markerCount) const
{
    if (markerCount == 0)
        return false;

    if (start == 0 && size == 0)
        return true;

    return start + size - 1 < markerCount;
}

MarkerIndexList MarkerSubset::toMarkerIndexList(unsigned int markerCount) const
{
    if (!isValid(markerCount))
        return {};

    const auto listSize = size == 0 ? markerCount : size;
    std::vector<unsigned int> markers(listSize);
    std::iota(markers.begin(), markers.end(), start);

    return markers;
}

bool isValid(const std::vector<MarkerSubset> subsets, unsigned int markerCount)
{
    if (subsets.empty()) {
        std::cout << "Subset list is empty" << std::endl;
        return false;
    }

    unsigned int sum = subsets.front().size;
    auto overlapCheck = [&sum](const MarkerSubset &a, const MarkerSubset &b) {
        sum += b.size;
        return a.last() != (b.first() - 1);
    };

    auto overlap = std::adjacent_find(subsets.cbegin(), subsets.cend(), overlapCheck);

    if (overlap != subsets.cend()) {
        const auto distance = std::distance(subsets.cbegin(), overlap);
        const auto next = overlap + 1;
        std::cout << "Subset " << distance << " lacks contiguity with the next subset:" << std::endl
                  << distance << ": [" << overlap->first() << ", " << overlap->last() << "]" << std::endl
                  << distance + 1 << ": [" << next->first() << ", " << next->last() << "]" << std::endl;
        return false;
    }

    if (sum != markerCount) {
        std::cout << "Expected " << markerCount << " markers, found: " << sum << std::endl;
        return false;
    }

    return true;
}

std::vector<MarkerSubset> generateEqualSubsets(unsigned int subsetCount, unsigned int markerCount)
{
    if (subsetCount == 0) {
        std::cout << "subsetCount cannot be 0" << std::endl;
        return {};
    }

    if (markerCount == 0) {
        std::cout << "markerCount cannot be 0" << std::endl;
        return {};
    }

    const auto size = static_cast<unsigned int>(std::ceil(static_cast<double>(markerCount) / static_cast<double>(subsetCount)));
    unsigned int start = 0;

    auto generateSubset = [&start, size]() {
        MarkerSubset s = {start, size};
        start = static_cast<unsigned int>(s.last()) + 1;
        return s;
    };

    std::vector<MarkerSubset> subsets(subsetCount);
    std::generate(subsets.begin(), subsets.end(), generateSubset);

    subsets.back().clamp(markerCount);

    return subsets;
}

bool operator==(const MarkerSubset &lhs, const MarkerSubset &rhs)
{
    return lhs.asTuple() == rhs.asTuple();
}

bool operator!=(const MarkerSubset &lhs, const MarkerSubset &rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& out, const MarkerSubset& subset)
{
    out.write(reinterpret_cast<const char *>(&subset), k_markerSubsetSize);
    return out;
}

std::istream &operator>>(std::istream &in, MarkerSubset &subset)
{
    in.read(reinterpret_cast<char *>(&subset), k_markerSubsetSize);
    return in;
}
