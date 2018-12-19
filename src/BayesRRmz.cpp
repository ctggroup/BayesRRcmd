/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "BayesRRmz.h"
#include "compression.h"
#include "data.hpp"
#include "distributions_boost.hpp"
#include "limitsequencegraph.h"
#include "options.hpp"
#include "parallelalgo.h"
#include "samplewriter.h"

#include <chrono>
#include <numeric>
#include <random>

BayesRRmz::BayesRRmz(Data &data, Options &opt)
    : flowGraph(new LimitSequenceGraph)
    , data(data)
    , opt(opt)
    , bedFile(opt.bedFile + ".bed")
    , outputFile(opt.mcmcSampleFile)
    , seed(opt.seed)
    , max_iterations(opt.chainLength)
    , burn_in(opt.burnin)
    , thinning(opt.thin)
    , dist(opt.seed)
    , usePreprocessedData(opt.analysisType == "PPBayes")
    , showDebug(false)
{
    float* ptr =static_cast<float*>(&opt.S[0]);
    cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();
}

BayesRRmz::~BayesRRmz()
{
}

void BayesRRmz::init(int K, unsigned int markerCount, unsigned int individualCount)
{
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
    mu = 0.0;       // mean or intercept
    sigmaG = 0.0;   // genetic variance
    sigmaE = 0.0;   // residuals variance

    // Linear model variables
    beta = VectorXd(markerCount);           // effect sizes
    y_tilde = VectorXd(individualCount);    // variable containing the adjusted residuals to exclude the effects of a given marker
    epsilon = VectorXd(individualCount);    // variable containing the residuals

    y = VectorXd();
    Cx = VectorXd();

    // Init the working variables
    const int km1 = K - 1;
    cVa[0] = 0;
    cVa.segment(1, km1) = cva;
    priorPi[0] = 0.5;
    priorPi.segment(1, km1) = priorPi[0] * cVa.segment(1, km1).array() / cVa.segment(1, km1).sum();
    y_tilde.setZero();

    cVaI[0] = 0;
    cVaI.segment(1, km1) = cVa.segment(1, km1).cwiseInverse();
    beta.setZero();
    sigmaG = dist.beta_rng(1,1);

    pi = priorPi;

    y = (data.y.cast<double>().array() - data.y.cast<double>().mean());
    y /= sqrt(y.squaredNorm() / (double(individualCount - 1)));

    epsilon = (y).array() - mu;
    sigmaE = epsilon.squaredNorm() / individualCount * 0.5;
}

int BayesRRmz::runGibbs()
{
    const unsigned int M(data.numIncdSnps);
    const unsigned int N(data.numKeptInds);
    const double NM1 = double(N - 1);
    const int K(int(cva.size()) + 1);
    const int km1 = K - 1;

    init(K, M, N);

    SampleWriter writer;
    writer.setFileName(outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.open();

    // Sampler variables
    VectorXd sample(2*M+4+N); // varible containg a sambple of all variables in the model, M marker effects, M component assigned to markers, sigmaE, sigmaG, mu, iteration number and Explained variance
    std::vector<unsigned int> markerI(M);
    std::iota(markerI.begin(), markerI.end(), 0);

    std::cout << "Running Gibbs sampling" << endl;
    const auto t1 = std::chrono::high_resolution_clock::now();

    // We can use a single fixed size buffer to decompress the preprocessed
    // data into as each column is the same size when uncompressed.
    unsigned char *decompressBuffer = nullptr;
    Map<VectorXf> Cx(nullptr, 1);
    const unsigned int colSize = N * sizeof(float);
    decompressBuffer = new unsigned char[colSize];
    new (&Cx) Map<VectorXf>(reinterpret_cast<float *>(decompressBuffer), N);

    // This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
    VectorXd components(M);
    components.setZero();
    for (unsigned int iteration = 0; iteration < max_iterations; iteration++) {
        // Output progress
        const auto startTime = std::chrono::high_resolution_clock::now();
        //if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
            std::cout << "iteration " << iteration << ": ";

        const double sigmaEpsilon = parallelStepAndSumEpsilon(epsilon, mu);
        parallelStepMuEpsilon(mu, epsilon, sigmaEpsilon, double(N), sigmaE, dist);

        std::random_shuffle(markerI.begin(), markerI.end());

        m0 = 0;
        v.setZero();

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution
        for (unsigned int j = 0; j < M; j++) {
            double acum = 0.0;
            const auto marker = markerI[j];

            // TODO: Can we improve things by decompressing a compressed mmap datafile?
            //const VectorXf Cx = getSnpData(marker);
            //const VectorXf &Cx = data.mappedZ.col(marker);
            extractData(reinterpret_cast<unsigned char *>(data.ppBedMap) + data.ppbedIndex[marker].pos,
                        static_cast<unsigned int>(data.ppbedIndex[marker].size),
                        decompressBuffer,
                        colSize);

            // Now y_tilde = Y-mu - X * beta + X.col(marker) * beta(marker)_old
            parallelUpdateYTilde(y_tilde, epsilon, Cx, beta(marker));

            // muk for the zeroth component=0
            muk[0] = 0.0;

            // We compute the denominator in the variance expression to save computations
            const double sigmaEOverSigmaG = sigmaE / sigmaG;
            denom = NM1 + sigmaEOverSigmaG * cVaI.segment(1, km1).array();

            // We compute the dot product to save computations
            const double num = parallelDotProduct(Cx, y_tilde);

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

            // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
            parallelUpdateEpsilon(epsilon, y_tilde, Cx, beta(marker));
        }

        m0 = int(M) - int(v[0]);
        sigmaG = dist.inv_scaled_chisq_rng(v0G + m0, (beta.col(0).squaredNorm() * m0 + v0G * s02G) / (v0G + m0));

        if (showDebug)
            printDebugInfo();

        const double epsilonSqNorm = parallelSquaredNorm(epsilon);
        sigmaE = dist.inv_scaled_chisq_rng(v0E + N, (epsilonSqNorm + v0E * s02E) / (v0E + N));
        pi = dist.dirichilet_rng(v.array() + 1.0);

        if (iteration >= burn_in && iteration % thinning == 0) {
            sample << iteration, mu, beta, sigmaE, sigmaG, components, epsilon;
            writer.write(sample);
        }

        const auto endTime = std::chrono::high_resolution_clock::now();
        const auto iterationDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        std::cout << iterationDuration / double(1000.0) << "s" << std::endl;
    }

    const auto t2 = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "duration: " << duration << "s" << std::endl;

    return 0;
}

void BayesRRmz::printDebugInfo() const
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
