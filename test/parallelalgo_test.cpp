/*
 * parallelalgo_test.cpp
 *
 *  Created on: 9 Jan 2019
 *      Author: admin
 */

#include<gtest/gtest.h>
#include "parallelalgo.h"
#include<Eigen/Eigen>
#include<distributions_boost.hpp>
//we test that the parallel squared norm returns the same value as the squared norm
const int N=300000;
const int samples=100;
//test that parallel squared norm is equivalent to vectorised squared norm
TEST(ParallelalgoTest,testparallelSquaredNorm){
	VectorXd input(N);
	input.setRandom();
	const double par=parallelSquaredNorm(input);
	const double vec=input.squaredNorm();
	ASSERT_NEAR (par, vec,1e-8);
}
//test that the parallel dot product is equivalent to the vectorised dot product
TEST(ParallelalgoTest,testparallelDotProduct){
	VectorXd y_tilde(N);
    VectorXd column(N);
	y_tilde.setRandom();
	column.setRandom();
	const double par=parallelDotProduct(column,y_tilde);
    const double vec=column.dot(y_tilde);
	ASSERT_NEAR (par, vec,1e-8);
}

//test thtat the parallel sum(epsilon-mu) is equivalent to the parallel verion
TEST(ParallelalgoTest,testparallelStepAndSumEpsilon){
	VectorXd epsilon(N);
	const double mu=0.1;
	epsilon.setRandom();

    const double vec=(epsilon.array()+mu).sum();
    const double par=parallelStepAndSumEpsilon(epsilon,mu);
    ASSERT_NEAR (par, vec,1e-8);
}
//test that mu is updated by the parallel funtion
TEST(ParallelalgoTest,testparallelStepMuEpsilon_muUpdated){
	VectorXd epsilon(N);
	double mu_old;
    double mu;
	const double sigmaE=0.5;
	double sigmaEpsilon;
	Distributions_boost dist(123);
	epsilon.setRandom();
	mu_old=0.1;
	mu=mu_old;
	sigmaEpsilon=(epsilon.array()-mu).sum();
	parallelStepMuEpsilon(mu, epsilon, sigmaEpsilon,  N, sigmaE,dist);
	ASSERT_NE (mu_old, mu);
}

//test that epsilon is updated by the parallel funtion
TEST(ParallelalgoTest,testparallelStepMuEpsilon_epsilonUpdated){
	VectorXd epsilon(N);
	double mu_old;
	double mu;
	VectorXd epsilon_expect(N);
	const double sigmaE=0.5;
	double sigmaEpsilon;
	Distributions_boost dist(123);
	Distributions_boost dist_2(123);
	epsilon.Random();
	epsilon_expect=epsilon;
	sigmaEpsilon=(epsilon.array()).sum();
	mu_old=dist.norm_rng(sigmaEpsilon/(double)N, sigmaE/(double)N);

	epsilon_expect=epsilon.array()-mu_old;
	mu=0.1;
	parallelStepMuEpsilon(mu, epsilon, sigmaEpsilon,  N, sigmaE,dist_2);
	ASSERT_TRUE (epsilon.isApprox(epsilon_expect));
}
//test that mu comes from the correct distribution.
TEST(ParallelalgoTest,testparallelStepMuEpsilon_muisNormal){
	VectorXd epsilon(N);
	VectorXd muVec(samples);
	 double mu_old;
	double mu;
	VectorXd epsilon_old(N);
	const double sigmaE=0.5;
	double sigmaEpsilon;
	double residualMean;
	Distributions_boost dist(123);
	epsilon_old.setRandom();
    mu_old=0.1;
    mu=mu_old;
	sigmaEpsilon=(epsilon.array()-mu_old).sum();
	for( int i=0;i<samples;i++){
		epsilon=epsilon_old;
		parallelStepMuEpsilon(mu, epsilon, sigmaEpsilon,  N, sigmaE,dist);
		muVec[i]=mu;
	}
	residualMean=muVec.mean()-(epsilon_old.sum()/(double) N);
    ASSERT_TRUE(residualMean<1e-6);
}

TEST(ParallelalgoTest,parallelUpdateYTilde){
	VectorXd epsilon(N);
    VectorXd column(N);
	VectorXd y_tilde(N);
	double beta;
	VectorXd y_tilde_expect(N);
	const double sigmaE=0.5;
	double sigmaEpsilon;
	Distributions_boost dist(123);
	epsilon.setRandom();
	column.setRandom();
	y_tilde.setRandom();
	beta=0.1;
    y_tilde_expect=epsilon+beta*column;
	parallelUpdateYTilde(y_tilde, epsilon,  column,  beta);
	ASSERT_TRUE (y_tilde.isApprox(y_tilde_expect));
}
TEST(ParallelalgoTest,parallelUpdateEpsilon){
	VectorXd epsilon(N);
    VectorXd column(N);
	VectorXd y_tilde(N);
	double beta;
	VectorXd epsilon_expect(N);
	const double sigmaE=0.5;
	double sigmaEpsilon;
	Distributions_boost dist(123);
	epsilon.setRandom();
	column.setRandom();
	y_tilde.setRandom();
	beta=0.1;
    epsilon_expect=y_tilde-beta*column;
	parallelUpdateEpsilon(epsilon, y_tilde,  column,  beta);
	ASSERT_TRUE (epsilon.isApprox(epsilon_expect));
}
