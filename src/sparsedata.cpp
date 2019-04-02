#include "sparsedata.h"

#include "compression.h"

SparseData::SparseData()
    : Data()
{

}

SparseData::~SparseData()
{

}

void SparseData::readBedFileSparse(const string &bedFile)
{
    cout << "Processing sparse data structures for bed file: " << bedFile << endl;
    if (numSnps == 0)
        throw ("Error: No SNP is retained for analysis.");
    if (numInds == 0)
        throw ("Error: No individual is retained for analysis.");

    ifstream BIT(bedFile.c_str(), ios::binary);
    if (!BIT)
        throw ("Error: can not open the file [" + bedFile + "] to read.");

    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    BIT.read(header, 3);
    if (!BIT || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01)
        throw ("Error: Incorrect first three bytes of bed file: " + bedFile);

    // Initialise data vectors with correct size
    means.setZero(numSnps);
    sds.setZero(numSnps);
    sqrdZ.setZero(numSnps);
    Zsum.setZero(numSnps);

    initialise();

    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 10: hetezygote; 01: missing
    for (unsigned int snp = 0; snp < numSnps; snp++) {
        SnpInfo *snpInfo = snpInfoVec[snp];

        const unsigned int size = (numInds + 3) >> 2;
        if (!snpInfo->included) {
            BIT.ignore(size);
            continue;
        }

        beginSnpColumn(snp);

        // Record the number of missing genotypes
        unsigned int missingGenotypeCount = 0;

        for (unsigned int i = 0; i < numInds;) {
            char ch;
            BIT.read(&ch, 1);
            if (!BIT)
                throw ("Error: problem with the BED file ... has the FAM/BIM file been changed?");

            bitset<8> b = ch;
            unsigned int k = 0;

            while (k < 7 && i < numInds) {
                if (!indInfoVec[i]->kept) {
                    k += 2;
                } else {
                    const unsigned int allele1 = (!b[k++]);
                    const unsigned int allele2 = (!b[k++]);

                    processAllele(snp, i, allele1, allele2);

                    if (allele1 == 0 && allele2 == 1) {  // missing genotype
                        ++missingGenotypeCount;
                    } else if (allele1 == 1 || allele2 == 1) { // Not zero
                        // Populate data for 1 or 2
                        const double value = allele1 + allele2;
                        means[snp] += value;
                        sqrdZ[snp] += value * value;
                        Zsum[snp] += value;
                    }
                }
                i++;
            }
        }

        // Calculate mean
        means[snp] /= static_cast<double>(numInds);

        // Call endSnpColumn before we calculate the standard deviation to allow us to
        // impute missing values if required.
        endSnpColumn(snp, missingGenotypeCount);

        // Calculate sds
        const double mean = means[snp];
        sds[snp] = std::sqrt((sqrdZ[snp] - 2 * mean * Zsum[snp] + static_cast<double>(numInds) * mean * mean) /
                             (static_cast<double>(numInds) - 1.0));

    }

    BIT.clear();
    BIT.close();

    cout << "Genotype data for " << numInds << " individuals and " << numSnps << " SNPs are included from [" + bedFile + "]." << endl;
}

double SparseData::computeNum(const unsigned int marker,
                              const double beta_old,
                              const VectorXd &epsilon,
                              const double epsilonSum) const
{
    return beta_old * (static_cast<double>(numInds) - 1.0) - means(marker) * epsilonSum / sds(marker) + dot(marker, epsilon);
}

double SparseData::computeEpsilonSumUpdate(const unsigned int marker,
                                           const double beta_old,
                                           const double beta) const
{
    //Regardless of which scheme, the update of epsilonSum is the same
    const double dBeta = beta_old - beta;
    return dBeta * Zsum(marker) / sds(marker) - dBeta * means(marker) * static_cast<double>(numInds) / sds(marker);
}

bool SparseData::writeStatistics(std::ofstream &outStream) const
{
    if (outStream.fail()) {
        std::cerr << "Error: unable to write SparseData statistics!" << std::endl;
        return false;
    }

    outStream.write(reinterpret_cast<const char *>(&numSnps), sizeof(numSnps));

    const std::streamsize size = sizeof(double) * numSnps;
    outStream.write(reinterpret_cast<const char *>(&means[0]), size);
    outStream.write(reinterpret_cast<const char *>(&sds[0]), size);
    outStream.write(reinterpret_cast<const char *>(&sqrdZ[0]), size);
    outStream.write(reinterpret_cast<const char *>(&Zsum[0]), size);
    outStream.flush();

    return true;
}

unsigned long SparseData::writeStatisticsCompressed(ofstream &outStream,
                                                    ofstream &indexStream) const
{
    if (outStream.fail()) {
        std::cerr << "Error: unable to write compressed SparseData statistics!" << std::endl;
        return 0;
    }

    if (indexStream.fail()) {
        std::cerr << "Error: unable to write compressed SparseData statistics index!" << std::endl;
        return 0;
    }

    unsigned long pos = sizeof(numSnps);
    outStream.write(reinterpret_cast<const char *>(&numSnps), pos);

    const auto maxCompressedOutputSize = maxCompressedDataSize<double>(numInds);
    unsigned char *compressedBuffer = new unsigned char[maxCompressedOutputSize];

    compressAndWriteWithIndex(means, outStream, indexStream, pos, compressedBuffer, maxCompressedOutputSize);
    compressAndWriteWithIndex(sds, outStream, indexStream, pos, compressedBuffer, maxCompressedOutputSize);
    compressAndWriteWithIndex(sqrdZ, outStream, indexStream, pos, compressedBuffer, maxCompressedOutputSize);
    compressAndWriteWithIndex(Zsum, outStream, indexStream, pos, compressedBuffer, maxCompressedOutputSize);

    delete[] compressedBuffer;
    return pos;
}
