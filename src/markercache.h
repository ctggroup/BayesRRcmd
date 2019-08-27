#ifndef MARKERCACHE_H
#define MARKERCACHE_H

#include "common.h"

#include <memory>
#include <vector>

class Data;
class Options;

class MarkerCache {
public:
    void clear();
    void populate(const Data *data, const Options *options);
    ConstMarkerPtr marker(unsigned int i) const;

protected:
    using MarkerPtrList = std::vector<ConstMarkerPtr>;
    MarkerPtrList m_markers;
};

const std::shared_ptr<MarkerCache> markerCache();

#endif // MARKERCACHE_H
