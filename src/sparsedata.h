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


    //OR
    //std::vector<std:vector<int>> Zones(M)//vector containing the vectors the indexes of elements of the bed matrix which are one for each column
    //std::vector<std:vector<int>> Ztwos(M)//vector containing the vectors the indexes of elements of the bed matrix which are two for each column

    //This following function, reads the bed file and stores each columns mean, sd, sum of squares and sum of elements, additionally it saves
    //in memory either the sparse matrix Zg or the two vector of vectors of indexes.
    void readBedFileSparse(const string &bedFile);

protected:
    virtual void initialise() = 0;
    virtual void beginSnpColumn(unsigned int snp) = 0;
    virtual void processAllele(unsigned int snp, unsigned int individual, unsigned int allele1, unsigned int allele2) = 0;
    virtual void endSnpColumn(unsigned int snp, unsigned int missingGenotypeCount) = 0;
};

#endif // SPARSEDATA_H
