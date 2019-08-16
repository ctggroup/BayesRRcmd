#ifndef MARKERBUILDER_H
#define MARKERBUILDER_H

#include "marker.h"

#include <memory>
#include <vector>

struct IndexEntry;

class MarkerBuilder
{
public:
    virtual ~MarkerBuilder();

    virtual void initialise(const unsigned int snp,
                            const unsigned int numInds);

    virtual void processAllele(unsigned int individual,
                               unsigned int allele1,
                               unsigned int allele2) = 0;

    virtual void endColumn() = 0;

    virtual void read(const std::string &file, const IndexEntry &index) const;

    virtual void decompress(unsigned char *data,
                            const IndexEntry &index) const;

    virtual std::unique_ptr<Marker> build();

protected:
    std::unique_ptr<Marker> m_marker = nullptr;

    unsigned int m_snp = 0;
    double m_numInds = 0;
    std::vector<unsigned int> m_missingIndices;
    double m_sum = 0;

    virtual void reset(const double numInds = 0);
    virtual void initialiseMarker();
};

#endif // MARKERBUILDER_H
