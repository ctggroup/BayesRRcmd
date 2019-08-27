#ifndef SPARSEMARKER_H
#define SPARSEMARKER_H

#include "marker.h"

struct SparseMarker : public Marker
{
    double mean = 0;
    double sd = 0;
    double sqrdZ= 0;
    double Zsum = 0;

    virtual void updateStatistics(unsigned int allele1, unsigned int allele2);

    std::streamsize size() const override;
    void read(std::istream *inStream) override;
    void write(std::ostream *outStream) const override;
};

#endif // SPARSEMARKER_H
