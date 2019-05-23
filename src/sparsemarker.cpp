#include "sparsemarker.h"

#include "compression.h"

#include <fstream>
#include <iostream>

double SparseMarker::computeNum(VectorXd &epsilon, const double beta_old)
{
    return computeNum(epsilon, beta_old, epsilonSum);
}

VectorXdPtr SparseMarker::calculateEpsilonChange(const double beta_old, const double beta)
{
    // now every update only saves delta epsilon sum
    epsilonSum = computeEpsilonSumUpdate(beta_old, beta);

    // return nullptr because we don't calculate the epsilon change here
    return nullptr;
}

void SparseMarker::updateStatistics(unsigned int allele1, unsigned int allele2)
{
    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        return;
    } else if (allele1 == 1 || allele2 == 1) { // Not zero
        // Populate data for 1 or 2
        const double value = allele1 + allele2;
        mean += value;
        sqrdZ += value * value;
        Zsum += value;
    }
}

std::streamsize SparseMarker::size() const
{
    return sizeof(double) * 4;
}

void SparseMarker::read(std::istream *inStream)
{
    if (inStream->fail()) {
        std::cerr << "Error: unable to read SparseMarker statistics!" << std::endl;
        return;
    }

    auto readDouble = [&](double *d) {
        inStream->read(reinterpret_cast<char *>(d), sizeof(double));
    };

    readDouble(&mean);
    readDouble(&sd);
    readDouble(&sqrdZ);
    readDouble(&Zsum);
}

void SparseMarker::write(std::ostream *outStream) const
{
    if (outStream->fail()) {
        std::cerr << "Error: unable to write SparseMarker statistics!" << std::endl;
        return;
    }

    auto writeDouble = [&](const double d) {
        outStream->write(reinterpret_cast<const char *>(&d), sizeof(double));
    };

    writeDouble(mean);
    writeDouble(sd);
    writeDouble(sqrdZ);
    writeDouble(Zsum);
}

double SparseMarker::computeNum(VectorXd &epsilon, const double beta_old, const double epsilonSum)
{
    return beta_old * (static_cast<double>(numInds) - 1.0) - mean * epsilonSum / sd + dot(epsilon);
}

double SparseMarker::computeEpsilonSumUpdate(const double beta_old, const double beta) const
{
    //Regardless of which scheme, the update of epsilonSum is the same
    const double dBeta = beta_old - beta;
    return dBeta * Zsum / sd - dBeta * mean * static_cast<double>(numInds) / sd;
}
