#include "raggedsparsedata.h"

#include <iterator>

#include "compression.h"

RaggedSparseData::RaggedSparseData()
    : SparseData()
{

}

double RaggedSparseData::dot(const unsigned int marker,
                             const VectorXd &epsilon) const
{
    return (epsilon(Zones[marker]).sum() + 2 * epsilon(Ztwos[marker]).sum()) / sds(marker);
}

void RaggedSparseData::updateEpsilon(VectorXd &epsilon,
                                     const unsigned int marker,
                                     const double beta_old,
                                     const double beta) const
{
    const double dBeta = beta_old - beta;
    const auto meanAdjustment = dBeta * means(marker) / sds(marker);
    // 1. Adjust for the means. If snp is 0, this will be the only adjustment made
    epsilon.array() -= meanAdjustment;

    // 2. Adjust for snp 1 values
    const double oneAdjustment = dBeta / sds(marker);
    epsilon(Zones[marker]).array() += oneAdjustment;

    // 3. Adjust for snp 2 values
    epsilon(Ztwos[marker]).array() += 2 * oneAdjustment;

    // 4. For missing values, undo step 1
    epsilon(Zmissing[marker]).array() += meanAdjustment;
}

bool RaggedSparseData::writeSparseData(const std::string &outFile,
                                       const bool compressed) const
{
    std::ofstream outStream(outFile.c_str(), std::ios::binary);
    if (outStream.fail()) {
        std::cerr << "Error: unable to open the output file for writing: " << outFile << std::endl;
        return false;
    }

    if (compressed) {
        const auto indexFile = outFile + "index";
        std::ofstream indexStream(indexFile.c_str(), std::ios::binary);
        if (indexStream.fail()) {
            std::cerr << "Error: unable to open the output index file for writing: " << indexFile << std::endl;
            return false;
        }

        unsigned long position = writeStatisticsCompressed(outStream, indexStream);
        if (position == 0)
            return false;

        return writeRaggedVectorCompressed(Zones, outStream, indexStream, position) &&
                writeRaggedVectorCompressed(Ztwos, outStream, indexStream, position) &&
                writeRaggedVectorCompressed(Zmissing, outStream, indexStream, position);
    } else {
        return writeStatistics(outStream) &&
                writeRaggedVector(Zones, outStream) &&
                writeRaggedVector(Ztwos, outStream) &&
                writeRaggedVector(Zmissing, outStream);
    }
}

void RaggedSparseData::initialise()
{
    // Initialise our RaggedVectors with IndexVectors with the estimated amount of reserved space.
    using size_type = RaggedVector::size_type;
    size_type estimatedDataCount = static_cast<size_type>(static_cast<double>(numInds) * 0.2);
    IndexVector tmp;
    tmp.reserve(estimatedDataCount);

    Zones = RaggedVector(numSnps, tmp);
    Ztwos = RaggedVector(numSnps, tmp);

    size_type estimatedMissingCount = static_cast<size_type>(static_cast<double>(numInds) * 0.01);
    tmp.reserve(estimatedMissingCount);

    Zmissing = RaggedVector(numSnps, tmp);
}

void RaggedSparseData::beginSnpColumn(unsigned int snp)
{
    m_currentOnes = &Zones[snp];
    m_currentTwos = &Ztwos[snp];
    m_currentMissing = &Zmissing[snp];
}

void RaggedSparseData::processAllele(unsigned int snp,
                                     unsigned int individual,
                                     unsigned int allele1,
                                     unsigned int allele2)
{
    assert(m_currentOnes);
    assert(m_currentTwos);
    assert(m_currentMissing);

    (void) snp; // Unused

    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        m_currentMissing->emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 0) {
        m_currentOnes->emplace_back(individual);
    } else if (allele1 == 1 && allele2 == 1) {
        m_currentTwos->emplace_back(individual);
    }
}

void RaggedSparseData::endSnpColumn(unsigned int snp,
                                    unsigned int missingGenotypeCount)
{
    (void) snp; // Unused
    (void) missingGenotypeCount; // Unused

    // No need to impute missing values

    m_currentOnes = nullptr;
    m_currentTwos = nullptr;
    m_currentMissing = nullptr;
}

bool RaggedSparseData::writeRaggedVector(const RaggedSparseData::RaggedVector &vector,
                                         ostream &outStream) const
{
    if (outStream.fail()) {
        std::cerr << "Error: unable to write RaggedVector!" << std::endl;
        return false;
    }

    // Assumes numSnps has been written in SparseData::writeStatistics
    for (const auto& indices : vector) {
        const IndexVector::size_type size = indices.size();
        outStream.write(reinterpret_cast<const char *>(&size), sizeof(IndexVector::size_type));

        if (size > 0)
            std::copy(indices.cbegin(), indices.cend(), std::ostream_iterator<IndexVector::value_type>(outStream));
    }

    outStream.flush();
    return true;
}

bool RaggedSparseData::writeRaggedVectorCompressed(const RaggedSparseData::RaggedVector &vector,
                                                   ostream &outStream,
                                                   ostream &indexStream,
                                                   unsigned long &position) const
{
    if (outStream.fail()) {
        std::cerr << "Error: unable to write compressed RaggedVector!" << std::endl;
        return false;
    }

    if (indexStream.fail()) {
        std::cerr << "Error: unable to write compressed RaggedVector index!" << std::endl;
        return false;
    }

    const auto maxCompressedOutputSize = maxCompressedDataSize<IndexVector::value_type>(numInds);
    unsigned char *compressedBuffer = new unsigned char[maxCompressedOutputSize];

    // Assumes numSnps has been written in SparseData::writeStatistics
    for (const auto& indices : vector) {
        compressAndWriteWithIndex<IndexVector::value_type>(indices,
                                                           outStream,
                                                           indexStream,
                                                           position,
                                                           compressedBuffer,
                                                           maxCompressedOutputSize);
    }

    delete[] compressedBuffer;
    return true;
}
