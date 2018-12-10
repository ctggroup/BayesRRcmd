/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "BayesRRm.h"
#include "data.hpp"
#include "distributions_boost.hpp"
#include "options.hpp"
#include "parallelalgo.h"
#include "samplewriter.h"

#include <chrono>
#include <numeric>
#include <random>
#include <Eigen/Eigen>

BayesRRm::BayesRRm(Data &data, Options &opt, const long memPageSize)
    : LinearReg(opt,data)
{
    float* ptr =static_cast<float*>(&opt.S[0]);
    cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();
    K=int(cva.size()) + 1;
    // Component variables
      priorPi = VectorXd(K);      // prior probabilities for each component
      pi = VectorXd(K);           // mixture probabilities
      cVa = VectorXd(K);          // component-specific variance
      logL = VectorXd(K);         // log likelihood of component
      muk = VectorXd (K);         // mean of k-th component marker effect size
      denom = VectorXd(K - 1);    // temporal variable for computing the inflation of the effect variance for a given non-zero componnet
      m0 = 0;                     // total num ber of markes in model
      v = VectorXd(K);            // variable storing the component assignment
      cVaI = VectorXd(K);         // inverse of the component variances

      // Mean and residual variables
      sigmaG = 0.0;   // genetic variance


      Cx = VectorXd();

      // Init the working variables
      cVa[0] = 0;
      cVa.segment(1, km1) = cva;
      priorPi[0] = 0.5;
      priorPi.segment(1, km1) = priorPi[0] * cVa.segment(1, km1).array() / cVa.segment(1, km1).sum();

      cVaI[0] = 0;
      cVaI.segment(1, km1) = cVa.segment(1, km1).cwiseInverse();
      sigmaG = dist.beta_rng(1,1);

      pi = priorPi;


      components=VectorXd(M);
      components.setZero();
      // sample vector and marker index vector
      sample=VectorXd(2*M+4+N);

      NM1 = double(N - 1);
      km1 = K - 1;


      m0 = 0;
      v.setZero();

      init_header();
}

BayesRRm::~BayesRRm()
{
}

void init()
{
	 // Sampler variables

}

void BayesRRm::updateBetas(int marker,const VectorXf &Cx)
{
            double acum = 0.0;


            // TODO: Can we improve things by decompressing a compressed mmap datafile?


            // muk for the zeroth component=0
            muk[0] = 0.0;

            // We compute the denominator in the variance expression to save computations
            const double sigmaEOverSigmaG = sigmaE / sigmaG;
            denom = NM1 + sigmaEOverSigmaG * cVaI.segment(1, km1).array();

            // We compute the dot product to save computations using Thanasis fix
            const double num = parallelcwiseProduct(Cx, y_tilde);

            // muk for the other components is computed according to equaitons
            muk.segment(1, km1) = num / denom.array();

            // Update the log likelihood for each component
            const double logLScale = sigmaG / sigmaE * NM1;
            logL = pi.array().log(); // First component probabilities remain unchanged
            logL.segment(1, km1) = logL.segment(1, km1).array()
                    - 0.5 * ((logLScale * cVa.segment(1, km1).array() + 1).array().log())
                    + 0.5 * (muk.segment(1, km1).array() * num) / sigmaE;

            double p(dist.beta_rng(1,1)); //I use beta(1,1) because I cant be bothered in using the std::random or create my own uniform distribution, I will change it later

            if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) {
                acum = 0;
            } else {
                acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
            }

            for (int k = 0; k < K; k++) {
                if (p <= acum) {
                    //if zeroth component
                    if (k == 0) {
                        beta(marker) = 0;
                    } else {
                        beta(marker) = dist.norm_rng(muk[k], sigmaE/denom[k-1]);
                    }
                    v[k] += 1.0;
                    components[marker] = k;
                    break;
                } else {
                    //if too big or too small
                    if (((logL.segment(1, km1).array() - logL[k+1]).abs().array() > 700).any()) {
                        acum += 0;
                    } else {
                        acum += 1.0 / ((logL.array() - logL[k+1]).exp().sum());
                    }
                }
            }
            betasqn+=beta(marker)*beta(marker); //cheaper to compute here than  again in the update hyper

}


void BayesRRm::updateHyper(){
	  m0 = int(M) - int(v[0]);
	// std::cout<<"beta norm "<<beta.squaredNorm()<<"\n";
	 //std::cout<<"v "<<v<<"\n";
	 //std::cout<< "v0g*s02g"<<v0G * s02G;
	 sigmaG = dist.inv_scaled_chisq_rng(v0G + m0, (betasqn * m0 + v0G * s02G) / (v0G + m0));
     std::cout<<"sigmaG "<<sigmaG << "\n";
	 pi = dist.dirichilet_rng(v.array() + 1.0);
     v.setZero();//restart component specific counters
     betasqn=0;
}

void BayesRRm::collectSample(){
	 sample << mu, beta, sigmaE, sigmaG, components, epsilon;
}

VectorXd BayesRRm::getSnpData(unsigned int marker) const
{
    if (!usePreprocessedData) {
        // Compute the SNP data length in bytes
        const size_t snpLenByt = (data.numInds % 4) ? data.numInds / 4 + 1 : data.numInds / 4;

        // I use a temporary variable to do the cast, there should be better ways to do this.
        VectorXf normedSnpData(data.numKeptInds);
        data.getSnpDataFromBedFileUsingMmap_openmp(bedFile, snpLenByt, memPageSize, marker, normedSnpData);
        return normedSnpData.cast<double>();
    } else {
        const VectorXf &sourceData = data.mappedZ.col(marker);
        VectorXd result(sourceData.size());
        parallelCastDouble(sourceData, result);
        return result;
    }
}

 void BayesRRm::printDebugInfo()
{
    const unsigned int N(data.numKeptInds);
    cout << "inv scaled parameters " << v0G + m0 << "__" << (beta.squaredNorm() * m0 + v0G * s02G) / (v0G + m0);
    cout << "num components: " << opt.S.size();
    cout << "\nMixture components: " << cva[0] << " " << cva[1] << " " << cva[2] << "\n";
    cout << "sigmaG: " << sigmaG << "\n";
    cout << "y mean: " << y.mean() << "\n";
    cout << "y sd: " << sqrt(y.squaredNorm() / (double(N - 1))) << "\n";
    cout << "x mean " << Cx.mean() << "\n";
    cout << "x sd " << sqrt(Cx.squaredNorm() / (double(N - 1))) << "\n";
}

 std::string& BayesRRm::getHeader(){
	return header;
}

 void BayesRRm::init_header(){
	 header="";
	 header+="mu,";
	     for (unsigned int i = 0; i < M; ++i) {
	         header+= "beta[" ;
	         header+=std::to_string(i+1);
	         header+= "],";
	     }

	     header+= "sigmaE,";
	     header+="sigmaG,";
	     for (unsigned int i = 0; i < M; ++i) {
	         header+="comp[";
	         header+=std::to_string(i+1);
	         header+="],";
	     }

	     unsigned int i;
	     for (i = 0; i < (N - 1); ++i) {
	         header += "epsilon[";
	         header += std::to_string(i + 1);
	         header+="],";
	     }

	     header+= "epsilon[";
	     header+= std::to_string(i + 1);
	     header+="]";
	     header+= "\n";
	//TODO header function
}


