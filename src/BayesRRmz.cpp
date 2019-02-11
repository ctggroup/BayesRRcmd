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
    float* ptr =static_cast<float*>(&opt.S[0]);
    m_cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();

    if (opt.analysisType == "PPAsyncBayes") {
        std::cerr << "PPAsyncBayes not implemented yet. Aborting" << std::endl;
        throw "PPAsyncBayes not implemented yet. Aborting";
        //flowGraph.reset(new AsyncSequenceGraph(this));
    } else {
        m_flowGraph.reset(new LimitSequenceGraph(this));
    }
}

BayesRRmz::~BayesRRmz()
{
}

void BayesRRmz::init(int K, unsigned int markerCount, unsigned int individualCount)
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

    m_y = (m_data.y.cast<double>().array() - m_data.y.cast<double>().mean());
    m_y /= sqrt(m_y.squaredNorm() / (double(individualCount - 1)));

    m_epsilon = (m_y).array() - m_mu;
    m_sigmaE = m_epsilon.squaredNorm() / individualCount * 0.5;
}

int BayesRRmz::runGibbs()
{
    const unsigned int M(m_data.numIncdSnps);
    const unsigned int N(m_data.numKeptInds);
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

    for (unsigned int iteration = 0; iteration < m_maxIterations; iteration++) {
        // Output progress
        const auto startTime = std::chrono::high_resolution_clock::now();
        //if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
            std::cout << "iteration " << iteration << ": ";

        //
#ifdef PARUP
            const double sigmaEpsilon = parallelStepAndSumEpsilon(m_epsilon, m_mu);
            parallelStepMuEpsilon(m_mu, m_epsilon, sigmaEpsilon, double(N), m_sigmaE, m_dist);
#else
         epsilon = epsilon.array() + mu;//  we substract previous value
         mu = dist.norm_rng(epsilon.sum() / (double)N, sigmaE / (double)N); //update mu
         epsilon = epsilon.array() - mu;// we substract again now epsilon =Y-mu-X*beta
#endif
        std::random_shuffle(markerI.begin(), markerI.end());

        m_m0 = 0;
        m_v.setZero();

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution.
        // The flow graph is constructed to allow the data to be decompressed in parallel for enforce sequential processing of each column
        // in turn. HOwever, within each column we make use of Intel TBB's parallel_for to parallelise the operations on the large vectors
        // of data.
        m_flowGraph->exec(N, M, markerI);

        m_m0 = int(M) - int(m_v[0]);
        m_sigmaG = m_dist.inv_scaled_chisq_rng(m_v0G + m_m0, (m_betasqn * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0));

        if (m_showDebug)
            printDebugInfo();
#ifdef PARUP
        const double epsilonSqNorm = parallelSquaredNorm(m_epsilon);
#else
        const double epsilonSqNorm = epsilon.squaredNorm();
#endif
        m_sigmaE = m_dist.inv_scaled_chisq_rng(m_v0E + N, (epsilonSqNorm + m_v0E * m_s02E) / (m_v0E + N));
        m_pi = m_dist.dirichilet_rng(m_v.array() + 1.0);

        if (iteration >= m_burnIn && iteration % m_thinning == 0) {
            sample << iteration, m_mu, m_beta, m_sigmaE, m_sigmaG, m_components, m_epsilon;
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

void BayesRRmz::processColumn(unsigned int marker, const Map<VectorXd> &Cx)
{
    const unsigned int N(m_data.numKeptInds);
    const double NM1 = double(N - 1);
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    double acum = 0.0;
    double beta_old;

    beta_old = m_beta(marker);

    // Now y_tilde = Y-mu - X * beta + X.col(marker) * beta(marker)_old
    if (m_components(marker) != 0.0) {
#ifdef PARUP
        parallelUpdateYTilde(m_y_tilde, m_epsilon, Cx, m_beta(marker));
#else
        y_tilde = epsilon + beta_old * Cx;
#endif
    } else {
        m_y_tilde = m_epsilon;
    }
    // muk for the zeroth component=0
    m_muk[0] = 0.0;

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;
    m_denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();

    // We compute the dot product to save computations
    // We compute the dot product to save computations
#ifdef PARUP
      const double num = parallelDotProduct(Cx, m_y_tilde);
#else
      const double num = Cx.dot(y_tilde);
#endif
    // muk for the other components is computed according to equaitons
    m_muk.segment(1, km1) = num / m_denom.array();

    // Update the log likelihood for each component
    VectorXd logL(K);
    const double logLScale = m_sigmaG / m_sigmaE * NM1;
    logL = m_pi.array().log(); // First component probabilities remain unchanged
    logL.segment(1, km1) = logL.segment(1, km1).array()
            - 0.5 * ((logLScale * m_cVa.segment(1, km1).array() + 1).array().log())
            + 0.5 * (m_muk.segment(1, km1).array() * num) / m_sigmaE;

    double p(m_dist.unif_rng());

    if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) {
        acum = 0;
    } else {
        acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
    }

    for (int k = 0; k < K; k++) {
        if (p <= acum) {
            //if zeroth component
            if (k == 0) {
                m_beta(marker) = 0;
            } else {
                m_beta(marker) = m_dist.norm_rng(m_muk[k], m_sigmaE/m_denom[k-1]);
            }
            m_v[k] += 1.0;
            m_components[marker] = k;
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
    m_betasqn += m_beta(marker) * m_beta(marker) - beta_old * beta_old;

    if (m_components(marker) != 0.0) {
#ifdef PARUP
        parallelUpdateEpsilon(m_epsilon, m_y_tilde, Cx, m_beta(marker));
#else
        epsilon = y_tilde - beta(marker) * Cx;
#endif
    } else {
        m_epsilon = m_y_tilde;
    }
    // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
}

void BayesRRmz::processColumnAsync(unsigned int marker, const Map<VectorXd> &Cx)
{
    // Lock and take local copies of needed variabls
    // [*] m_beta(marker) rwr - used, updated, then used - per column, could take a copy and update at end
    // [*] m_betasqn w - updated here, used in BayezRRmz::runGibbs
    // [*] m_components(marker) rwr - used, updated, then used - per column, could take a copy and update at end
    // [ ] m_y_tilde wr - updated first, then used throughout
    // [ ] m_epsilon rw - used throughout, then updated
    // [ ] m_v w - updated here, used in BayezRRmz::runGibbs

    // Temporaries
    // - cost of locking vs allocating per iteration?
    // m_denom wr - computed from m_cVaI
    // m_muk wr - computed from m_cVaI

    // m_data.numKeptInds r - could be a member?
    // m_cva.size() r - could be a member?
    // m_sigmaE r - calculated in BayesRRmz::init
    // m_sigmaG r - calculated in BayesRRmz::init, updated in BayezRRmz::runGibbs
    // m_pi r - calculated in BayesRRmz::init, updated in BayezRRmz::runGibbs
    // m_dist r - the engine is not thread safe
    // m_cVa r - calculated in BayesRRmz::init
    // m_cVaI r - calculated in BayesRRmz::init

    double beta = 0;
    double component = 0;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        beta = m_beta(marker);
        component = m_components(marker);
    }
    const double beta_old = beta;

    // Do work

    // Now y_tilde = Y-mu - X * beta + X.col(marker) * beta(marker)_old
    if (component != 0.0) {
#ifdef PARUP
        parallelUpdateYTilde(m_y_tilde, m_epsilon, Cx, beta);
#else
        y_tilde = epsilon + beta_old * Cx;
#endif
    } else {
        m_y_tilde = m_epsilon;
    }
    // muk for the zeroth component=0
    m_muk[0] = 0.0;

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;

    const double NM1 = double(m_data.numKeptInds - 1);
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    m_denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();

    // We compute the dot product to save computations
    // We compute the dot product to save computations
#ifdef PARUP
      const double num = parallelDotProduct(Cx, m_y_tilde);
#else
      const double num = Cx.dot(y_tilde);
#endif
    // muk for the other components is computed according to equaitons
    m_muk.segment(1, km1) = num / m_denom.array();

    // Update the log likelihood for each component
    VectorXd logL(K);
    const double logLScale = m_sigmaG / m_sigmaE * NM1;
    logL = m_pi.array().log(); // First component probabilities remain unchanged
    logL.segment(1, km1) = logL.segment(1, km1).array()
            - 0.5 * ((logLScale * m_cVa.segment(1, km1).array() + 1).array().log())
            + 0.5 * (m_muk.segment(1, km1).array() * num) / m_sigmaE;


    double acum = 0.0;
    if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) {
        acum = 0;
    } else {
        acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
    }

    const double p(m_dist.unif_rng());
    for (int k = 0; k < K; k++) {
        if (p <= acum) {
            //if zeroth component
            if (k == 0) {
                beta = 0;
            } else {
                beta = m_dist.norm_rng(m_muk[k], m_sigmaE/m_denom[k-1]);
            }
            m_v[k] += 1.0;
            component = k;
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

    if (component != 0.0) {
#ifdef PARUP
        parallelUpdateEpsilon(m_epsilon, m_y_tilde, Cx, beta);
#else
        epsilon = y_tilde - beta(marker) * Cx;
#endif
    } else {
        m_epsilon = m_y_tilde;
    }
    // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new

    // Lock to write updates (at end, or perhaps as updates are computed)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_beta(marker) = beta;
        m_betasqn += beta * beta - beta_old * beta_old;
        m_components(marker) = component;
    }
}

void BayesRRmz::printDebugInfo() const
{
    const unsigned int N(m_data.numKeptInds);
    cout << "inv scaled parameters " << m_v0G + m_m0 << "__" << (m_beta.squaredNorm() * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0);
    cout << "num components: " << m_opt.S.size();
    cout << "\nMixture components: " << m_cva[0] << " " << m_cva[1] << " " << m_cva[2] << "\n";
    cout << "sigmaG: " << m_sigmaG << "\n";
    cout << "y mean: " << m_y.mean() << "\n";
    cout << "y sd: " << sqrt(m_y.squaredNorm() / (double(N - 1))) << "\n";
//    cout << "x mean " << Cx.mean() << "\n";
//    cout << "x sd " << sqrt(Cx.squaredNorm() / (double(N - 1))) << "\n";
}
