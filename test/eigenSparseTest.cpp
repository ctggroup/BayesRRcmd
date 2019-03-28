#include<gtest/gtest.h>
#include<Eigen/Eigen>
#include <distributions_boost.hpp>
#include "data.hpp"
#include "options.hpp"
#include "eigensparsedata.h"



class EigenSparseDataTest public ::testing::Test {
 protected:
  void SetUp() override{
   data = DataPtr(new EigenSparseData);
   data->readFamFile(famFile);
   data->readBimFile(bimFile);
   data->readPhenotypeFile(bedFile);
   epsilon.resize();
   y_tilde.resize();
   epsilon.setRandom();
   y_tilde.setRandom();
   beta=0.1;
   beta_old=0.2;
}
  string bedFile;
  string phenotypeFile;
  string famFile;
  string bimFile;
  DataPtr data;
  const VectorXd epsilon;
  const VectorXd y_tilde;
  const double beta;
  const double beta_old;
}

/*function to test

void EigenSparseData::updateEpsilon(VectorXd &epsilon, const unsigned int marker, const double beta_old, const double beta) const
{
    const double dBeta = beta_old - beta;
    epsilon += dBeta * Zg[marker] / sds(marker) - dBeta * means(marker) / sds(marker) * m_ones;
}

 */
TEST(EigenSDTest,testUpdateEpsilon){
  //lets test 10 markers
  //for(int  marker=0;marker<10; marker++){
  int marker=1;
  VectoXd Centered(data.Zg[marker] - data.Zg[marker].mean() );
  sds=Centered.squaredNorm()/(data->-1.0);	
  Centered/=sds;	
  updateExpected= epsilon+(beta_old-beta)*Centered;
  updateGot= data->updateEpsilon(epsilon, marker, beta_old,beta);
  EXPECT_TRUE(updateExpected.isAprox(updateGot));
    //}
}

