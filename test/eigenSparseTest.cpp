#include<gtest/gtest.h>
#include<Eigen/Eigen>
#include <distributions_boost.hpp>
#include "data.hpp"
#include "sparsedata.h"
#include "options.hpp"
#include "eigensparsedata.h"
typedef Eigen::Triplet<double> T;
using namespace Eigen;
/*Ideally we would test the code of eigensparsedata functions, but for some reason I cant call it
 * without producing a segfalt
 *
 *
 *
 */ 

class EigenSparseDataTests:public::testing::Test {
 protected:
 // string bedFile="../../BayesRRcmd/test/data/uk10k_chr1_1mb.bed";
 // string phenotypeFile="../../BayesRRcmd/test/data/test.phen";
 // string famFile="../../BayesRRcmd/test/data/uk10k_chr1_1mb.fam";
 // string bimFile="../../BayesRRcmd/test/data/uk10k_chr1_1mb.bim";
  const int numInds=100000;  
  EigenSparseData data;
  VectorXd epsilon;
  const double beta=0.1;
  const double beta_old=0.2;
  std::vector<SparseVector<double>> Zg;//sparse column
  double dBeta;
  VectorXd Z;//dense column
  double means;
  double sds;
  VectorXd Centered;
  double var;
  void SetUp() override{
  //    data.readFamFile(famFile);
  //    data.readBimFile(bimFile);
  //    data.readBedFileSparse(bedFile);
  //    data=EigenSparseData();
  //    data.numInds=100000;
      std::default_random_engine gen;
      std::uniform_real_distribution<double> dist(0.0,1.0);
      std::vector<T> tripletList;
      SparseVector<double> column(numInds);
      Z=VectorXd(numInds);
      Z.setZero();
     
      for(int i=0; i< numInds; i++)
      { 
             auto v_i=dist(gen);
             if(v_i<0.2){ 
                column.coeffRef(i)=(double)v_i;
                Z[i]=(double)v_i;
             }
      }
      data.Zg.push_back(column);      
      epsilon=VectorXd(data.numInds); 
      epsilon.setRandom();
      Zg.push_back(column);
      dBeta=beta_old-beta;
      means=column.sum()/((double)numInds);
      sds= sqrt((column.squaredNorm() - 2 * means * column.sum() + static_cast<double>(numInds) * means * means)/((double)numInds-1.0));     
      Centered=Z.array()-Z.mean();
      var=Centered.squaredNorm()/((double)numInds-1.0);
      Centered.array()/=sqrt(var);
   }
};

TEST_F(EigenSparseDataTests,testUpdateEpsilon){
     int marker=0;
     VectorXd updateExpected;
     VectorXd updateGot;
     VectorXd ones(numInds);
     ones.setOnes();
     double norm;
     updateExpected= (beta_old-beta)*Centered;
     
     updateGot=VectorXd(dBeta * Zg[marker] / sds)- dBeta * means / sds * ones;
     norm=(updateExpected-updateGot).squaredNorm();
     EXPECT_DOUBLE_EQ(norm,0.0);
}
    
   
