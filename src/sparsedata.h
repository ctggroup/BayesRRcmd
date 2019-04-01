#ifndef SPARSEDATA_H
#define SPARSEDATA_H

#include "data.hpp"

class SparseData : public Data
{
public:
    SparseData();
    virtual ~SparseData();

    VectorXd means; //vector that contains the mean of each column of the bed file matrix
    VectorXd sds; //vector that contains the sd of each column of the bed file matrix
    VectorXd sqrdZ; //vector that contains the sum of squares of each column of the bed file matrix
    VectorXd Zsum;  //vector that contains the sum of elements of each column of the bed file matrix

    void readBedFileSparse(const string &bedFile);

    virtual bool writeSparseData(const string &outFile, const bool compress) const = 0;

    virtual void updateEpsilon(VectorXd &epsilon, const unsigned int marker, const double beta_old, const double beta) const = 0;

    virtual double computeNum(const unsigned int marker, const double beta_old, const VectorXd &epsilon, const double epsilonSum) const;
    virtual double computeEpsilonSumUpdate(const unsigned int marker, const double beta_old, const double beta) const;

protected:
    virtual double dot(const unsigned int marker, const VectorXd &epsilon) const = 0;

    virtual void initialise() = 0;
    virtual void beginSnpColumn(unsigned int snp) = 0;
    virtual void processAllele(unsigned int snp, unsigned int individual, unsigned int allele1, unsigned int allele2) = 0;
    virtual void endSnpColumn(unsigned int snp, unsigned int missingGenotypeCount) = 0;

    bool writeStatistics(std::ofstream& outStream) const;
    unsigned long writeStatisticsCompressed(std::ofstream& outStream, ofstream &indexStream) const;
};

#endif // SPARSEDATA_H
