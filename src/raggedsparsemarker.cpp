#include "raggedsparsemarker.h"

#include <iostream>

std::streamsize RaggedSparseMarker::size() const
{
    const auto valueTypeSize = sizeof(RaggedSparseMarker::IndexVector::value_type);
    const auto sizeTypeSize = sizeof(RaggedSparseMarker::IndexVector::size_type);

    return SparseMarker::size() +
            static_cast<std::streamsize>(
                sizeTypeSize + (valueTypeSize * Zones.size()) +
                sizeTypeSize + (valueTypeSize * Ztwos.size()) +
                sizeTypeSize + (valueTypeSize * Zmissing.size()));
}

void RaggedSparseMarker::read(std::istream *inStream)
{
    if (inStream->fail()) {
        std::cerr << "Error: unable to read RaggedSparseMarker!" << std::endl;
        return;
    }

    SparseMarker::read(inStream);

    auto readIndexVector = [&](IndexVector &v) {
        IndexVector::size_type size = 0;
        inStream->read(reinterpret_cast<char *>(&size),
                         sizeof(IndexVector::size_type));

        v.clear();
        v.resize(size);
        if (size > 0)
            inStream->read(reinterpret_cast<char *>(v.data()),
                           static_cast<std::streamsize>(size * sizeof (IndexVector::value_type)));
    };

    readIndexVector(Zones);
    readIndexVector(Ztwos);
    readIndexVector(Zmissing);
}

void RaggedSparseMarker::write(std::ostream *outStream) const
{
    if (outStream->fail()) {
        std::cerr << "Error: unable to write RaggedSparseMarker!" << std::endl;
        return;
    }

    SparseMarker::write(outStream);

    using IndexVector = RaggedSparseMarker::IndexVector;
    auto writeIndexVector = [&](const IndexVector & v) {
        const IndexVector::size_type size = v.size();
        outStream->write(reinterpret_cast<const char *>(&size),
                         sizeof(IndexVector::size_type));

        if (size > 0)
            outStream->write(reinterpret_cast<const char *>(v.data()),
                             static_cast<std::streamsize>(size * sizeof (IndexVector::value_type)));
    };

    writeIndexVector(Zones);
    writeIndexVector(Ztwos);
    writeIndexVector(Zmissing);
}

bool RaggedSparseMarker::isValid() const
{
    const bool invalid = Zones.empty() && Ztwos.empty();
    if (invalid)
        std::cerr << "SNPs that do not vary are should be removed prior to analysis. "
                  << "Otherwise, this message indicates a decompression error"
                  << std::endl;

    return !invalid;
}
