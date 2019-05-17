/*
 * BayesRRg.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "BayesRRg.h"
#include "data.hpp"
#include "distributions_boost.hpp"
#include "options.hpp"
#include "samplewriter.h"

#include <chrono>
#include <numeric>
#include <random>

BayesRRg::BayesRRg(Data &data, Options &opt, const long memPageSize)
: data(data)
, opt(opt)
, bedFile(opt.bedFile + ".bed")
, memPageSize(memPageSize)
, outputFile(opt.mcmcSampleFile)
, seed(opt.seed)
, max_iterations(opt.chainLength)
, burn_in(opt.burnin)
, thinning(opt.thin)
, dist(opt.seed)
, usePreprocessedData(opt.analysisType == "PPBayes")
, showDebug(false)
{
//  float* ptr =static_cast<float*>(&opt.S[0]);
//  cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();
}

BayesRRg::~BayesRRg()
{
}

void BayesRRg::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    // Group component variables
    priorPi = MatrixXd(numberGroups,K); // prior probabilities for each component
    pi = MatrixXd(numberGroups,K);      // mixture probabilities
    cVa = VectorXd(K);          		// component-specific variance
    cVaI = VectorXd(K);         		// inverse of the component variances
    logL = VectorXd(K);         		// log likelihood of component
    muk = VectorXd (K);         		// mean of k-th component marker effect size
    denom = VectorXd(K - 1);   			// temporal variable for computing the inflation of the effect variance for a given non-zero component
    m0 = 0;                     		// total number of markers in model
    v = MatrixXd(numberGroups,K);		// variable storing the component assignment

    // Mean and residual variables
    mu = 0.0;       					// mean or intercept
    sigmaG = 0.0;   					// genetic variance
    sigmaE = 0.0;   					// residuals variance
    sigmaGG = VectorXd(numberGroups);	// group genetic variance

    // Linear model variables
    beta = VectorXd(markerCount);           // effect sizes
    y_tilde = VectorXd(individualCount);    // variable containing the adjusted residuals to exclude the effects of a given marker
    epsilon = VectorXd(individualCount);    // variable containing the residuals

    //phenotype vector
    y = VectorXd();

    //SNP column vector
    Cx = VectorXd();

    //number of groups vector
    betaAcum = VectorXd(numberGroups);

    //init the working variables: cva and cvaI
    const int km1 = K - 1; // TODO isn't this already defined before calling init?

    cVa[0] = 0;
    cVa.segment(1, km1) = cva.row(data.G(marker));

    cVaI[0] = 0;
    cVaI.segment(1, km1) = cVa.segment(1, km1).cwiseInverse();

    //vector with component class for each marker
    components=VectorXd(markerCount);
    components.setZero();

    //set priors for pi parameters
    for(int i=0; i < numberGroups; i++){
    	priorPi.row(i)(0)=0.5;
    	for(int k=1;k<K;k++){
    		priorPi.row(i)(k)=0.5/K;
    	}
    }

    //set all elements in y_tilde and beta to zero
    y_tilde.setZero();
    beta.setZero();

    //sample sigmaGG from beta distribution
   	for(int i=0; i<numberGroups;i++)
   		sigmaGG[i]=dist.beta_rng(1,1);

   	//assign value to pi
    pi = priorPi;

    //scale phenotype vector stored in data.y
    y = (data.y.cast<double>().array() - data.y.cast<double>().mean());
    y /= sqrt(y.squaredNorm() / (double(individualCount - 1)));

    //initialize epsilon vector as the phenotype vector
    epsilon = (y).array() - mu;
    sigmaE = epsilon.squaredNorm() / individualCount * 0.5;
    betasqn=0;
    epsilonsum=0;
    ytildesum=0;
}

int BayesRRg::runGibbs()
{

	// initialize number of marker, individuals, groups and their mixtures
    const unsigned int M(data.numSnps);
    const unsigned int N(data.numInds);
    const double NM1 = double(N - 1);
    const int K(data.mS.cols() + 1);
    const int km1 = K - 1;
    const unsigned int numberGroups(data.numGroups);
    MatrixXd cva(data.mS);

    //initialize variables with init member function
    init(K, M, N);

    //specify how to write samples
    SampleWriter writer;
    writer.setFileName(outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.openGroup(); //instead of writer.open()

    // Sampler variables
    VectorXd sample(2*M+3+numberGroups+N);
    // VectorXd sample(2*M+4+N); // variable containing a sample of all variables in the model, M marker effects, M component assigned to markers, sigmaE, sigmaG, mu, iteration number and Explained variance
    std::vector<unsigned int> markerI(M);
    std::iota(markerI.begin(), markerI.end(), 0);

    srd::cout << "Number of groups: " << numberGroups << endl;

    std::cout << "Running Gibbs sampling" << endl;
    const auto t1 = std::chrono::high_resolution_clock::now();

    for (unsigned int iteration = 0; iteration < max_iterations; iteration++) {
        // Output progress
        const auto iterStart = std::chrono::high_resolution_clock::now();
        if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
            std::cout << "iteration: " << iteration << std::endl;

        epsilon = epsilon.array() + mu;//  we substract previous value
        mu = dist.norm_rng(epsilon.sum() / (double)N, sigmaE / (double)N); //update mu
        epsilon = epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta

        std::random_shuffle(markerI.begin(), markerI.end());

        m0 = 0;
        v.setZero();
        betaAcum.setZero();

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
        for (unsigned int j = 0; j < M; j++) {

            double acum = 0.0;
            const auto marker = markerI[j];
            double beta_old = betaAcum(data.G(marker)); //changed beta(marker) by betaAcum(data.G(marker))

            //get sigmaG of marker group index
            sigmaG = sigmaGG[data.G(marker)];

            //read data for column with member function getSnpData
            Cx = getSnpData(marker);

            // residual update only if marker was previously included in model
            if(components(marker)!=0){
                y_tilde=epsilon+beta_old*Cx;
            }
            else{
                y_tilde=epsilon;
            }
            // muk for the zeroth component=0
            muk[0] = 0.0;

            // We compute the denominator in the variance expression to save computations
            const double sigmaEOverSigmaG = sigmaE / sigmaG;
            denom = NM1 + sigmaEOverSigmaG * cVaI.segment(1, km1).array();

            // We compute the dot product to save computations
            const double num = (Cx.dot(y_tilde));

            // muk for the other components is computed according to equations
            muk.segment(1, km1) = num / denom.array();

            // Update the log likelihood for each variance component
            const double logLScale = sigmaG / sigmaE * NM1;
            logL = pi.row(data.G(marker)).array().log(); // First component probabilities remain unchanged
            logL.segment(1, km1) = logL.segment(1, km1).array()
                    		        - 0.5 * ((logLScale * cVa.segment(1, km1).array() + 1).array().log())
                    		        + 0.5 * (muk.segment(1, km1).array() * num) / sigmaE;

            double p(dist.unif_rng());

            if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) { // this quantity will be exponentiated and go to huge values and we set to infinity, thus the reciprocal is set to zero
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
                        betaAcum(data.G(marker))+= pow(beta(marker),2);
                    }
                    v.row(data.G(marker))(k)+=1.0;
                    components[marker] = k;
                    break;
                } else {
                    if (((logL.segment(1, km1).array() - logL[k+1]).abs().array() > 700).any()) {//again controlling for large values for the exponent
                        acum += 0;
                    } else {
                        acum += 1.0 / ((logL.array() - logL[k+1]).exp().sum());
                    }
                }
            }
         //   betasqn+=betaAcum(data.G(marker))*betaAcum(data.G(marker))-beta_old*beta_old; //changed beta(marker) by betaAcum(data.G(marker))
            // residual update only if updated marker is included in model
            // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
            if(components(marker)!=0){
                epsilon=y_tilde-betaAcum(data.G(marker))*Cx;
            }
            else{
                epsilon=y_tilde;
            }
        }
        //set no. of markers included in the model
        m0 = int(M) - int(v[0]);

        const double epsilonSqNorm=epsilon.squaredNorm();
        //sample residual variance sigmaE from inverse gamma
        sigmaE = dist.inv_scaled_chisq_rng(v0E + N, (epsilonSqNorm + v0E * s02E) / (v0E + N));

        //sample sigmaG from inverse gamma and hyperparameter pi from dirichlet
        //for each group
        for(int i = 0; i < numberGroups; i++){
        	m0 = v.row(i).sum() - v.row(i)(0);
        	sigmaGG[i] = dist.inv_scaled_chisq_rng(v0G + m0, (betaAcum(i) * m0 + v0G * s02G) / (v0G + m0));
        	// sigmaG = dist.inv_scaled_chisq_rng(v0G + m0, (betasqn * m0 + v0G * s02G) / (v0G + m0))
        	pi.row(i) = dist.dirichilet_rng(v.row(i).array() + 1.0);
        }


        if (showDebug)
            printDebugInfo();

        //write samples
        if (iteration >= burn_in && iteration % thinning == 0) {
            sample << iteration, mu, beta, sigmaE, sigmaGG, components, epsilon;
            writer.write(sample);
        }

        //output time taken for each iteration
        const auto endTime = std::chrono::high_resolution_clock::now();
        const auto dif = endTime - iterStart;
        const auto iterationDuration = std::chrono::duration_cast<std::chrono::milliseconds>(dif).count();
        std::cout << iterationDuration / double(1000.0) << "s" << std::endl;

        //end of iteration
    }
    //show info on total time spent
    const auto t2 = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "duration: " << duration << "s" << std::endl;

    return 0;
}

VectorXd BayesRRg::getSnpData(unsigned int marker) const
{
    if (!usePreprocessedData) {
        //read column from RAM loaded genotype matrix.
        return data.Z.col(marker).cast<double>();
    } else {
        //read column from preprocessed and memory mapped genotype matrix file.
        return data.mappedZ.col(marker).cast<double>();
    }
}

void BayesRRg::printDebugInfo() const
{
    const unsigned int N(data.numInds);
    cout << "inv scaled parameters " << v0G + m0 << "__" << (beta.squaredNorm() * m0 + v0G * s02G) / (v0G + m0);
    cout << "num components: " << data.mS.cols();
    cout << "\nMixture components:\n" << cva << "\n";
    cout << "sigmaG\n: " << sigmaGG << "\n";
    cout << "y mean: " << y.mean() << "\n";
    cout << "y sd: " << sqrt(y.squaredNorm() / (double(N - 1))) << "\n";
    // cout << "x mean " << Cx.mean() << "\n";
    // cout << "x sd " << sqrt(Cx.squaredNorm() / (double(N - 1))) << "\n";
}
