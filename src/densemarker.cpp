#include "densemarker.h"

#include "common.h"
#include "compression.h"

#include <fstream>

CompressedMarker DenseMarker::compress() const
{
    CompressedMarker compressed;
    const auto maxCompressedOutputSize = maxCompressedDataSize<double>(numInds);
    compressed.buffer.reset(new unsigned char[maxCompressedOutputSize]);
    compressed.index.originalSize = static_cast<unsigned long>(size());
    compressed.index.compressedSize = compressData(*Cx,
                                                   compressed.buffer.get(),
                                                   maxCompressedOutputSize);
    return compressed;
}

void DenseMarker::decompress(unsigned char *data, const IndexEntry &index)
{
    const unsigned int colSize = numInds * sizeof(double);
    buffer.reset(new unsigned char[colSize]);

    extractData(data + index.pos,
                static_cast<unsigned int>(index.compressedSize),
                buffer.get(),
                colSize);

    Cx = std::make_shared<Map<VectorXd>>(reinterpret_cast<double *>(buffer.get()),
                                         numInds);
}

std::streamsize DenseMarker::size() const
{
    return numInds * sizeof(double);
}

void DenseMarker::read(std::istream *inStream)
{
    inStream->read(reinterpret_cast<char *>(buffer.get()),
                   size());
}

void DenseMarker::write(std::ostream *outStream) const
{
    outStream->write(reinterpret_cast<char *>(buffer.get()),
                     size());
}

bool DenseMarker::isValid() const
{
    return buffer && Cx;
}
