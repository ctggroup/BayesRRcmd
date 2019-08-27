#ifndef RAGGEDSPARSEMARKER_H
#define RAGGEDSPARSEMARKER_H

#include "sparsemarker.h"
#include "markerbuilder.h"

struct RaggedSparseMarker : public SparseMarker
{
    using IndexVector = std::vector<int>;
    // the indexes of elements of the bed matrix which are one for this column
    IndexVector Zones;
    // the indexes of elements of the bed matrix which are two for this column
    IndexVector Ztwos;

    // the indexes of elements of the bed matrix which are missing for this column
    IndexVector Zmissing;

    std::streamsize size() const override;
    void read(std::istream *inStream) override;
    void write(std::ostream *outStream) const override;

    bool isValid() const override;
};

#endif // RAGGEDSPARSEMARKER_H
