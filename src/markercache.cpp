#include "markercache.h"

#include "data.hpp"
#include "marker.h"
#include "markerbuilder.h"
#include "options.hpp"

const std::shared_ptr<MarkerCache> markerCache()
{
    static const auto cache = std::make_shared<MarkerCache>();
    return cache;
}

void MarkerCache::clear()
{
    m_markers.clear();
}

void MarkerCache::populate(const Data *data, const Options *options)
{
    assert(data);
    assert(options);

    if (!data || !options)
        return;

    clear();
    m_markers.resize(data->numSnps);

    std::unique_ptr<MarkerBuilder> builder{builderForType(options->preprocessDataType)};

    unsigned int snp = data->markerSubset().first();
    auto beginItr = m_markers.begin();
    std::advance(beginItr, snp);
    auto endItr = m_markers.begin();
    std::advance(endItr, data->markerSubset().last() + 1);
    std::for_each(beginItr, endItr, [&snp, &builder, &data, &options](ConstMarkerPtr &marker) {
        builder->initialise(snp, data->numInds);
        const auto index = data->ppbedIndex[snp];
        if (options->compress) {
            builder->decompress(reinterpret_cast<unsigned char*>(data->ppBedMap), index);
        } else {
            builder->read(ppFileForType(*options), index);
        }
        marker = builder->build();
        ++snp;
    });
    std::cout << "Cached " << snp << " markers" << std::endl;
}

ConstMarkerPtr MarkerCache::marker(unsigned int i) const
{
    if (i > m_markers.size()) {
        std::cerr << "Requesting marker out of bounds: " << i << ". "
                  << m_markers.size() << " cached markers." << std::endl;
        assert(false);
        return {};
    }

    return m_markers.at(i);
}
