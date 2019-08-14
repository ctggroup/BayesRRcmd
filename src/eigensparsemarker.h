#ifndef EIGENSPARSEMARKER_H
#define EIGENSPARSEMARKER_H

#include "sparsemarker.h"

struct EigenSparseMarker : public SparseMarker
{
    using UnitDataType = double;
    SparseVector<UnitDataType> Zg;

    std::streamsize size() const override;
    void read(std::istream *inStream) override;
    void write(std::ostream *outStream) const override;

    bool isValid() const override;
};

#endif // EIGENSPARSEMARKER_H
