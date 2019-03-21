/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "BayesRBase.hpp"
#include "samplewriter.h"
#include "analysisgraph.hpp"

#include <chrono>

BayesRBase::BayesRBase(const Data *data, Options &opt)
    : m_flowGraph(nullptr)
    , m_data(data)
    , m_opt(opt)
    , m_bedFile(opt.bedFile + ".bed")
    , m_outputFile(opt.mcmcSampleFile)
    , m_seed(opt.seed)
    , m_maxIterations(opt.chainLength)
    , m_burnIn(opt.burnin)
    , m_thinning(opt.thin)
    , m_dist(opt.seed)
    , m_usePreprocessedData(opt.analysisType == "PPBayes")
    , m_showDebug(false)
{
    assert(m_data);

    float* ptr =static_cast<float*>(&opt.S[0]);
    m_cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();
}

BayesRBase::~BayesRBase()
{
}

void BayesRBase::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    // Component variables
    m_priorPi = VectorXd(K);      // prior probabilities for each component
    m_pi = VectorXd(K);           // mixture probabilities
    m_cVa = VectorXd(K);          // component-specific variance
    m_muk = VectorXd (K);         // mean of k-th component marker effect size
    m_denom = VectorXd(K - 1);    // temporal variable for computing the inflation of the effect variance for a given non-zero componnet
    m_m0 = 0;                     // total num ber of markes in model
    m_v = VectorXd(K);            // variable storing the component assignment
    m_cVaI = VectorXd(K);         // inverse of the component variances

    // Mean and residual variables
    m_mu = 0.0;       // mean or intercept
    m_sigmaG = 0.0;   // genetic variance
    m_sigmaE = 0.0;   // residuals variance

    // Linear model variables
    m_beta = VectorXd(markerCount);           // effect sizes
    m_y_tilde = VectorXd(individualCount);    // variable containing the adjusted residuals to exclude the effects of a given marker
    m_epsilon = VectorXd(individualCount);    // variable containing the residuals
    m_async_epsilon = VectorXd(individualCount);

    m_y = VectorXd();
    //Cx = VectorXd();

    // Init the working variables
    const int km1 = K - 1;
    m_cVa[0] = 0;
    m_cVa.segment(1, km1) = m_cva;
    m_priorPi[0] = 0.5;
    m_priorPi.segment(1, km1) = m_priorPi[0] * m_cVa.segment(1, km1).array() / m_cVa.segment(1, km1).sum();
    m_y_tilde.setZero();

    m_cVaI[0] = 0;
    m_cVaI.segment(1, km1) = m_cVa.segment(1, km1).cwiseInverse();
    m_beta.setZero();
    m_sigmaG = m_dist.beta_rng(1,1);

    m_pi = m_priorPi;

    m_y = (m_data->y.cast<double>().array() - m_data->y.cast<double>().mean());
    m_y /= sqrt(m_y.squaredNorm() / (double(individualCount - 1)));

    m_epsilon = (m_y).array() - m_mu;
    m_sigmaE = m_epsilon.squaredNorm() / individualCount * 0.5;
}

int BayesRBase::runGibbs()
{
    if (!m_flowGraph) {
        std::cout << "Cannot run Gibbs analysis without a flow graph!" << std::endl;
        return 1;
    }

    const unsigned int M(m_data->numSnps);
    const unsigned int N(m_data->numInds);
    const int K(int(m_cva.size()) + 1);

    init(K, M, N);

    SampleWriter writer;
    writer.setFileName(m_outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.open();

    // Sampler variables
    VectorXd sample(2*M+4+N); // varible containg a sambple of all variables in the model, M marker effects, M component assigned to markers, sigmaE, sigmaG, mu, iteration number and Explained variance
    std::vector<unsigned int> markerI(M);
    std::iota(markerI.begin(), markerI.end(), 0);

    std::cout << "Running Gibbs sampling" << endl;
    const auto t1 = std::chrono::high_resolution_clock::now();

    // This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
    m_components.resize(M);
    m_components.setZero();

    long meanIterationTime = 0;
    long meanFlowGraphIterationTime = 0;

    for (unsigned int iteration = 0; iteration < m_maxIterations; iteration++) {
        // Output progress
        const auto startTime = std::chrono::high_resolution_clock::now();
        //if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
        std::cout << "iteration " << iteration << ": ";

        m_epsilon = m_epsilon.array() + m_mu;//  we substract previous value
        m_mu = m_dist.norm_rng(m_epsilon.sum() / (double)N, m_sigmaE / (double)N); //update mu
        m_epsilon = m_epsilon.array() - m_mu;// we substract again now epsilon =Y-mu-X*beta

        if (m_isAsync)
            std::memcpy(m_async_epsilon.data(), m_epsilon.data(), static_cast<size_t>(m_epsilon.size()) * sizeof(double));

        std::random_shuffle(markerI.begin(), markerI.end());

        m_m0 = 0;
        m_v.setZero();

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution.
        // The flow graph is constructed to allow the data to be decompressed in parallel for enforce sequential processing of each column
        // in turn. HOwever, within each column we make use of Intel TBB's parallel_for to parallelise the operations on the large vectors
        // of data.
        const auto flowGraphStartTime = std::chrono::high_resolution_clock::now();
        m_flowGraph->exec(N, M, markerI);
        const auto flowGraphEndTime = std::chrono::high_resolution_clock::now();

        m_m0 = int(M) - int(m_v[0]);
        m_sigmaG = m_dist.inv_scaled_chisq_rng(m_v0G + m_m0, (m_betasqn * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0));

        if (m_showDebug)
            printDebugInfo();
        const double epsilonSqNorm = m_epsilon.squaredNorm();
        m_sigmaE = m_dist.inv_scaled_chisq_rng(m_v0E + N, (epsilonSqNorm + m_v0E * m_s02E) / (m_v0E + N));
        m_pi = m_dist.dirichilet_rng(m_v.array() + 1.0);

        if (iteration >= m_burnIn && iteration % m_thinning == 0) {
            sample << iteration, m_mu, m_beta, m_sigmaE, m_sigmaG, m_components, m_epsilon;
            writer.write(sample);
        }

        const auto endTime = std::chrono::high_resolution_clock::now();
        const auto iterationDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        const auto flowGraphDuration = std::chrono::duration_cast<std::chrono::milliseconds>(flowGraphEndTime - flowGraphStartTime).count();
        std::cout << static_cast<double>(iterationDuration) / 1000.0 << "s (" << static_cast<double>(flowGraphDuration) / 1000.0 << "s)" << std::endl;
        meanIterationTime += iterationDuration;
        meanFlowGraphIterationTime += flowGraphDuration;
    }

    const auto t2 = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "duration: " << duration << "s" << std::endl;
    const double meanIterationDuration = (static_cast<double>(meanIterationTime) / 1000.0) / static_cast<double>(m_maxIterations);
    const double meanFlowGraphIterationDuration = (static_cast<double>(meanFlowGraphIterationTime) / 1000.0) / static_cast<double>(m_maxIterations);
    std::cout << "mean iteration duration: " << meanIterationDuration  << "s" << std::endl
              << "mean flowgraph duration: " << meanFlowGraphIterationDuration << "s" << std::endl;

    return 0;
}

void BayesRBase::printDebugInfo() const
{
    const unsigned int N(m_data->numInds);
    cout << "inv scaled parameters " << m_v0G + m_m0 << "__" << (m_beta.squaredNorm() * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0);
    cout << "num components: " << m_opt.S.size();
    cout << "\nMixture components: " << m_cva[0] << " " << m_cva[1] << " " << m_cva[2] << "\n";
    cout << "sigmaG: " << m_sigmaG << "\n";
    cout << "y mean: " << m_y.mean() << "\n";
    cout << "y sd: " << sqrt(m_y.squaredNorm() / (double(N - 1))) << "\n";
    //    cout << "x mean " << Cx.mean() << "\n";
    //    cout << "x sd " << sqrt(Cx.squaredNorm() / (double(N - 1))) << "\n";
}
