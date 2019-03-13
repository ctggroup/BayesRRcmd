#include "data.hpp"

class DataG:Data{
    VectorXd means; //vector that contains the mean of each column of the bed file matrix
    VectorXd sds; //vector that contains the sd of each column of the bed file matrix
    VectorXd sqrdZ; //vector that contains the sum of squares of each column the bed file matrix
    VectorXd Zsum;  //vector that contains the sum of squares of each columnof the bed file matrix
    //SparseMatrix<double> Zg; //this is a sparse matrix that contains the uncentered and unscaled elements of the bed matrix
    //OR
    //std::vector<std:vector<int>> Zones(M)//vector containing the vectors the indexes of elements of the bed matrix which are one for each column
    //std::vector<std:vector<int>> Ztwos(M)//vector containing the vectors the indexes of elements of the bed matrix which are two for each column
    //
    //This following function, reads the bed file and stores each columns mean, sd, sum of squares and sum of elements, additionally it saves
    //in memory either the sparse matrix Zg or the two vector of vectors of indexes.
    void readBedFile_G(const string &bedFile);
}
