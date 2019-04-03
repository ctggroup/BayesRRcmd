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

double RaggedSparseMarker::dot(const VectorXd &epsilon) const
{
    return (epsilon(Zones).sum() + 2 * epsilon(Ztwos).sum()) / sd;
}

template<>
void MarkerBuilder<RaggedSparseMarker>::initialise(const unsigned int snp,
                                                   const unsigned int numInds)
{
    reset(numInds);

    m_marker.reset(new RaggedSparseMarker);
    m_marker->i = snp;
    m_marker->numInds = numInds;

    using size_type = RaggedSparseMarker::IndexVector::size_type;
    const size_type estimatedDataCount = static_cast<size_type>(numInds * 0.2);
    m_marker->Zones.reserve(estimatedDataCount);
    m_marker->Ztwos.reserve(estimatedDataCount);

    const size_type estimatedMissingCount = static_cast<size_type>(numInds * 0.01);
    m_marker->Zmissing.reserve(estimatedMissingCount);
}

template<>
void MarkerBuilder<RaggedSparseMarker>::processAllele(unsigned int individual,
                                                      unsigned int allele1,
                                                      unsigned int allele2)
{
    updateStatistics(m_marker.get(), allele1, allele2);

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        m_marker->Zmissing.emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 0) {
        m_marker->Zones.emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 1) {
        m_marker->Ztwos.emplace_back(individual);
    }
}

template<>
void MarkerBuilder<RaggedSparseMarker>::endColumn()
{
    // Calculate mean
    m_marker->mean /= static_cast<double>(m_marker->Zmissing.size());

    // Calculate sd
    const double mean = m_marker->mean;
    m_marker->sd = std::sqrt((m_marker->sqrdZ - 2.0 * mean * m_marker->Zsum + m_numInds * mean * mean) /
                         (m_numInds - 1.0));
}


template<>
CompressedMarker compress(const RaggedSparseMarker *marker)
{
    // Prepare a stream to write into
    const auto valueTypeSize = sizeof(RaggedSparseMarker::IndexVector::value_type);
    const auto sizeTypeSize = sizeof(RaggedSparseMarker::IndexVector::size_type);

    const auto bufferSize = statisticsSize() +
            sizeTypeSize + (valueTypeSize * marker->Zones.size()) +
            sizeTypeSize + (valueTypeSize * marker->Ztwos.size()) +
            sizeTypeSize + (valueTypeSize * marker->Zmissing.size());

    std::unique_ptr<char[]> buffer;
    buffer.reset(new char[bufferSize]);

    std::ostringstream stream;
    stream.rdbuf()->pubsetbuf(buffer.get(), bufferSize);

    // Write the marker to the stream
    write(marker, &stream);

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

template<>
void write(const RaggedSparseMarker *marker, std::ostream *outStream)
{
    if (!marker)
        return;

    if (outStream->fail()) {
        std::cerr << "Error: unable to write RaggedSparseMarker!" << std::endl;
        return;
    }

    if (!writeStatistics(marker, outStream))
        return;

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

    writeIndexVector(marker->Zones);
    writeIndexVector(marker->Ztwos);
    writeIndexVector(marker->Zmissing);
}
