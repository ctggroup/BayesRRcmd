#include "densemarker.h"

#include "common.h"
#include "compression.h"

#include <fstream>

double DenseMarker::computeNum(VectorXd &epsilon, const double beta_old)
{
    //in order to not break async and sync updates for dense we change this
    //we now CX dot CX = N-1 given that both are already centered and scaled
    return Cx->dot(epsilon) + beta_old * static_cast<double>(numInds-1);
}

VectorXdPtr DenseMarker::calculateEpsilonChange(const double beta_old, const double beta)
{
    return std::make_unique<VectorXd>((beta_old-beta) * *Cx);
}

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
