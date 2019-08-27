#ifndef data_hpp
#define data_hpp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <set>
#include <bitset>
#include <iomanip>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <boost/format.hpp>
#include "gadgets.hpp"
#include "common.h"


using namespace std;
using namespace Eigen;

class SnpInfo {
public:
    const string ID;
    const string a1; // the referece allele
    const string a2; // the coded allele
    const int chrom;
    const float genPos;
    const int physPos;

    int index;
    int window;
    int windStart;  // for window surrounding the SNP
    int windSize;   // for window surrounding the SNP
    float af;       // allele frequency
    bool included;  // flag for inclusion in panel
    bool isQTL;     // for simulation

    SnpInfo(const int idx, const string &id, const string &allele1, const string &allele2,
            const int chr, const float gpos, const int ppos)
    : ID(id), index(idx), a1(allele1), a2(allele2), chrom(chr), genPos(gpos), physPos(ppos) {
        window = 0;
        windStart = -1;
        windSize  = 0;
        af = -1;
        included = true;
        isQTL = false;
    };

    void resetWindow(void) {windStart = -1; windSize = 0;};
};

class IndInfo {
public:
    const string famID;
    const string indID;
    const string catID;    // catenated family and individual ID
    const string fatherID;
    const string motherID;
    const int famFileOrder; // original fam file order
    const int sex;  // 1: male, 2: female

    int index;
    bool kept;

    float phenotype;

    VectorXf covariates;  // covariates for fixed effects

    IndInfo(const int idx, const string &fid, const string &pid, const string &dad, const string &mom, const int sex)
    : famID(fid), indID(pid), catID(fid+":"+pid), fatherID(dad), motherID(mom), index(idx), famFileOrder(idx), sex(sex) {
        phenotype = -9;
        kept = true;
    }
};

using PpBedIndex = std::vector<IndexEntry>;

class Data {
public:
    Data();

    // mmap related data
    int ppBedFd;
    double *ppBedMap;
    Map<MatrixXd> mappedZ;
    PpBedIndex ppbedIndex;

    // Original data
    MatrixXf X;              // coefficient matrix for fixed effects
    MatrixXf Z;              // coefficient matrix for SNP effects
    VectorXf D;              // 2pqn
    VectorXf y;              // phenotypes
    vector<int> G; 			 // groups
    VectorXd fail;			 // Failure indicator

    // Vectors for the sparse format solution
    std::vector<std::vector<int>> Zones; // Vector for each SNP: per SNP all the indices with 1 allele are written down
    std::vector<std::vector<int>> Ztwos; // Vector for each SNP: per SNP all the indices with 2 alleles are written down
    VectorXd means; //Mean for each SNP
    VectorXd sds;
    VectorXd mean_sd_ratio;


    MatrixXf XPX;            // X'X the MME lhs
    MatrixXf ZPX;            // Z'X the covariance matrix of SNPs and fixed effects
    VectorXf XPXdiag;        // X'X diagonal
    VectorXf ZPZdiag;        // Z'Z diagonal
    VectorXf XPy;            // X'y the MME rhs for fixed effects
    VectorXf ZPy;            // Z'y the MME rhs for snp effects

    VectorXf snp2pq;         // 2pq of SNPs
    VectorXf se;             // se from GWAS summary data
    VectorXf tss;            // total ss (ypy) for every SNP
    VectorXf b;              // beta from GWAS summary data
    VectorXf n;              // sample size for each SNP in GWAS

    vector<SnpInfo*> snpInfoVec;
    vector<IndInfo*> indInfoVec;

    map<string, SnpInfo*> snpInfoMap;
    map<string, IndInfo*> indInfoMap;

    vector<SnpInfo*> incdSnpInfoVec;
    vector<IndInfo*> keptIndInfoVec;

    vector<string> fixedEffectNames;
    vector<string> snpEffectNames;


    vector<bool> fullSnpFlag;
    vector<vector<SnpInfo*> > mldmVec;

    unsigned numFixedEffects = 0;
    unsigned numSnps;
    unsigned numInds;

    unsigned numGroups = 1; // number of groups

    void preprocessCSVFile(const string &csvFile, const string &preprocessedCSVFile, const string &preprovessedCSVIndexFile, bool compress);
    void mapPreprocessBedFile(const string &preprocessedBedFile);
    void unmapPreprocessedBedFile();

    void mapCompressedPreprocessBedFile(const string &preprocessedBedFile, const string &indexFile);
    void unmapCompressedPreprocessedBedFile();

    void readCSV(const string &filename, int cols);

    void readFamFile(const string &famFile);
    void readBimFile(const string &bimFile);
    void readBedFile_noMPI(const string &bedFile);
    void readBedFile_noMPI_unstandardised(const string &bedFile);

    void readPhenotypeFile(const string &phenFile);
    void readGroupFile(const string &groupFile);
    void readCSVFile(const string &csvFile);
    void readCSVPhenFile( const string &csvFile);
    void readmSFile(const string& mSfile);
    //BayesW variables
    void readFailureFile(const string &failureFile);
};

#endif /* data_hpp */
