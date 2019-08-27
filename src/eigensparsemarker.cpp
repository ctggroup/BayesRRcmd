#include "eigensparsemarker.h"

#include <iostream>

std::streamsize EigenSparseMarker::size() const
{
    const auto count = static_cast<unsigned long>(Zg.nonZeros());
    return SparseMarker::size() +
            static_cast<std::streamsize>(
                sizeof (SparseVector<UnitDataType>::Index) +
                (sizeof(SparseVector<UnitDataType>::Scalar) * count) +
                (sizeof(SparseVector<UnitDataType>::StorageIndex) * count));
}

void EigenSparseMarker::read(std::istream *inStream)
{
    if (inStream->fail()) {
        std::cerr << "Error: unable to read EigenSparseMarker!" << std::endl;
        return;
    }

    SparseMarker::read(inStream);

    Eigen::Index count = 0;
    inStream->read(reinterpret_cast<char *>(&count),
                   sizeof(Eigen::Index));

    Zg.resize(numInds); // Number of rows;
    Zg.resizeNonZeros(count);
    if (count > 0) {
        inStream->read(reinterpret_cast<char *>(Zg.valuePtr()),
                       count * static_cast<std::streamsize>(sizeof(SparseVector<UnitDataType>::Scalar)));
        inStream->read(reinterpret_cast<char *>(Zg.innerIndexPtr()),
                       count * static_cast<std::streamsize>(sizeof(SparseVector<UnitDataType>::StorageIndex)));
    }
}

void EigenSparseMarker::write(std::ostream *outStream) const
{
    if (outStream->fail()) {
        std::cerr << "Error: unable to write EigenSparseMarker!" << std::endl;
        return;
    }

    SparseMarker::write(outStream);

    const Eigen::Index count = Zg.nonZeros();
    outStream->write(reinterpret_cast<const char *>(&count),
                     sizeof(Eigen::Index));

    if (count > 0) {
        outStream->write(reinterpret_cast<const char *>(Zg.valuePtr()),
                         count * static_cast<std::streamsize>(sizeof(SparseVector<UnitDataType>::Scalar)));
        outStream->write(reinterpret_cast<const char *>(Zg.innerIndexPtr()),
                         count * static_cast<std::streamsize>(sizeof(SparseVector<UnitDataType>::StorageIndex)));
    }
}

bool EigenSparseMarker::isValid() const
{
    if (Zg.nonZeros() <= 0) {
        std::cerr << "SNPs that do not vary are should be removed prior to analysis. "
                  << "Otherwise, this message indicates a decompression error"
                  << std::endl;
        return false;
    }

    // nonZeros can be > 0 but size can be 0 if the SparseVector hasn't been
    // built correctly.
    return Zg.size() > 0;
}
