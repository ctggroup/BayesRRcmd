#include "raggedsparsemarker.h"

#include "compression.h"

#include <fstream>
#include <iostream>
#include <iterator>

void RaggedSparseMarker::updateEpsilon(VectorXd &epsilon, const double beta_old, const double beta)
{
    SparseMarker::updateEpsilon(epsilon, beta_old, beta);

    const double dBeta = beta_old - beta;
    const auto meanAdjustment = dBeta * mean / sd;
    // 1. Adjust for the means. If snp is 0, this will be the only adjustment made
    epsilon.array() -= meanAdjustment;

    // 2. Adjust for snp 1 values
    const double oneAdjustment = dBeta / sd;
    epsilon(Zones).array() += oneAdjustment;

    // 3. Adjust for snp 2 values
    epsilon(Ztwos).array() += 2 * oneAdjustment;

    // 4. For missing values, undo step 1
    epsilon(Zmissing).array() += meanAdjustment;
}

CompressedMarker RaggedSparseMarker::compress() const
{
    // Prepare a stream to write into
    const auto bufferSize = size();
    std::unique_ptr<char[]> buffer;
    buffer.reset(new char[bufferSize]);

    std::ostringstream stream;
    stream.rdbuf()->pubsetbuf(buffer.get(), bufferSize);

    // Write the marker to the stream
    write(&stream);

    // Compress the stream
    const auto maxCompressedOutputSize = maxCompressedDataSize<char>(bufferSize);

    CompressedMarker compressed;
    compressed.buffer.reset(new unsigned char[maxCompressedOutputSize]);
    compressed.size = compressData(buffer.get(),
                                   bufferSize,
                                   compressed.buffer.get(),
                                   maxCompressedOutputSize);
    return compressed;
}

size_t RaggedSparseMarker::size() const
{
    const auto valueTypeSize = sizeof(RaggedSparseMarker::IndexVector::value_type);
    const auto sizeTypeSize = sizeof(RaggedSparseMarker::IndexVector::size_type);

    return SparseMarker::size() +
            sizeTypeSize + (valueTypeSize * Zones.size()) +
            sizeTypeSize + (valueTypeSize * Ztwos.size()) +
            sizeTypeSize + (valueTypeSize * Zmissing.size());
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
            std::copy(v.cbegin(),
                      v.cend(),
                      std::ostream_iterator<IndexVector::value_type>(*outStream));
    };

    writeIndexVector(Zones);
    writeIndexVector(Ztwos);
    writeIndexVector(Zmissing);
}

double RaggedSparseMarker::dot(const VectorXd &epsilon) const
{
    return (epsilon(Zones).sum() + 2 * epsilon(Ztwos).sum()) / sd;
}
