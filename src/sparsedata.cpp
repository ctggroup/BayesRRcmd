#include "sparsedata.h"

using T = std::tuple<Eigen::Index, SparseData::UnitDataType>;
using TupleList = std::vector<T>;

const SparseData::UnitDataType SparseData::kMissing = std::numeric_limits<SparseData::UnitDataType>::lowest();

SparseData::SparseData()
    : Data()
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

    Zg.resize(numSnps);

    // Data is expected to be about 80% zeros, so reserve a bit more than 20% of the expected space
    const TupleList::size_type estimatedDataCount = static_cast<TupleList::size_type>(static_cast<double>(numInds) * 0.25);

    // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 10: hetezygote; 01: missing
    for (unsigned int snp = 0; snp < numSnps; snp++) {
        SnpInfo *snpInfo = snpInfoVec[snp];

        const unsigned int size = (numInds + 3) >> 2;
        if (!snpInfo->included) {
            BIT.ignore(size);
            continue;
        }

        // Create the triplet list for our sparse data representation
        TupleList tuples;
        tuples.reserve(estimatedDataCount);

        // Record the number of missing genotypes
        auto missingGenotypeCount = 0;

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

                    if (allele1 == 0 && allele2 == 1) {  // missing genotype
                        tuples.emplace_back(i, kMissing);
                        ++missingGenotypeCount;
                    } else if (allele1 == 1 || allele2 == 1) { // Not zero
                        // Populate data for 1 or 2
                        const double value = allele1 + allele2;
                        tuples.emplace_back(i, static_cast<UnitDataType>(value));
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
        const double mean = means[snp];

        // Update sqrdZ with imputed values
        sqrdZ[snp] += (mean * mean) * missingGenotypeCount;

        // Update Zsum with imputed values
        Zsum[snp] += mean * missingGenotypeCount;

        // Create the SparseVector
        auto &vector = Zg.at(snp);
        vector.resize(numInds); // Number of rows
        vector.reserve(tuples.size()); // Number of rows that are not zero

        // Create a vector of doubles for sds computation
        VectorXd doubles;
        doubles.resize(numInds);
        doubles.fill(0 - mean);

        // Fill the SparseVector and doubles
        std::for_each(tuples.cbegin(), tuples.cend(), [&](const T &t) {
            const auto index = std::get<0>(t);
            auto value = std::get<1>(t);
            // If the value is zero or smaller, impute using the mean.
            if (value <= 0)
                value = mean;

            vector.insertBack(index) = value;
            doubles(index) = value - mean;
        });

        // Calculate sds
        sds[snp] = std::sqrt(doubles.squaredNorm() / (static_cast<double>(numInds) - 1.0));
    }

    BIT.clear();
    BIT.close();

    cout << "Genotype data for " << numInds << " individuals and " << numSnps << " SNPs are included from [" + bedFile + "]." << endl;
}
