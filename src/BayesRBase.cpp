/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "BayesRBase.hpp"
#include "bayesrkernel.h"
#include "samplewriter.h"
#include "analysisgraph.hpp"
#include "marker.h"
#include "logwriter.h"

#include <chrono>
#include <mutex>

BayesRBase::BayesRBase(const Data *data, const Options *opt)
    : Analysis(data, opt)
    , m_outputFile(opt->mcmcSampleFile)
    , m_iterLogFile(opt->iterLogFile)
    , m_seed(opt->seed)
    , m_maxIterations(opt->chainLength)
    , m_burnIn(opt->burnin)
    , m_thinning(opt->thin)
    , m_cva(opt->S)
    , m_dist(opt->seed)
    , m_showDebug(opt->iterLog)
    , m_colLog(opt->colLog)
    , m_colLogFile(opt->colLogFile)
    , m_v0E(opt->v0E)
    , m_s02E(opt->s02E)
    , m_v0G(opt->v0G)
    , m_s02G(opt->s02G)
{
    assert(m_data);
}

void BayesRBase::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    const auto groupCount = m_data->numGroups;
    // Component variables
    m_priorPi = MatrixXd(groupCount, K);    // prior probabilities for each component
    const auto priorPi = 0.5 / K;
    for (auto i = 0; i < m_priorPi.rows(); ++i) {
        m_priorPi.row(i).fill(priorPi);
        m_priorPi.row(i)(0) = 0.5;
    }

    m_pi = m_priorPi;                       // mixture probabilities
    m_cVa = VectorXd(K);                    // component-specific variance
    m_muk = VectorXd (K);                   // mean of k-th component marker effect size
    m_denom = VectorXd(K - 1);              // temporal variable for computing the inflation of the effect variance for a given non-zero componnet
    m_m0 = 0;                               // total num ber of markes in model
    m_v = MatrixXd(groupCount, K);          // variable storing the component assignment
    m_cVaI = VectorXd(K);                   // inverse of the component variances

    // Mean and residual variables
    m_mu = 0.0;       // mean or intercept
    m_sigmaG.resize(groupCount);
    std::generate(m_sigmaG.begin(), m_sigmaG.end(), [&dist = m_dist]() {
        return dist.unif_rng();
    });

    m_sigmaE = 0.0;   // residuals variance

    // Linear model variables
    m_beta = VectorXd::Zero(markerCount);           // effect sizes
    m_y_tilde = VectorXd::Zero(individualCount);    // variable containing the adjusted residuals to exclude the effects of a given marker
    m_epsilon = VectorXd(individualCount);    // variable containing the residuals
    m_betasqnG = VectorXd(groupCount);
    //acum values
    m_acum = VectorXd::Zero(markerCount);           // acum vector
    m_y = VectorXd();
    //Cx = VectorXd();

    // Init the working variables
    const int km1 = K - 1;
    m_cVa[0] = 0;
    m_cVa.segment(1, km1) = m_cva;

    m_cVaI[0] = 0;
    m_cVaI.segment(1, km1) = m_cVa.segment(1, km1).cwiseInverse();

    m_y = (m_data->y.cast<double>().array() - m_data->y.cast<double>().mean());
    m_y /= sqrt(m_y.squaredNorm() / (double(individualCount - 1)));

    m_epsilon = (m_y).array() - m_mu;
    m_sigmaE = m_epsilon.squaredNorm() / individualCount * 0.5;
    m_epsilonSum=m_epsilon.sum();

    m_randomNumbers.resize(markerCount);

    if(m_colLog)
    {
        m_colWriter.setFileName(m_colLogFile);
        m_colWriter.open();
    }
}

void BayesRBase::prepareForAnylsis()
{
    // Generate the random numbers required for this iteration. The random
    // number engine is not thread safe, so generate them up front to avoid
    // having to use a mutex.
    std::generate(m_randomNumbers.begin(), m_randomNumbers.end(), [&dist = m_dist]()
                  -> std::array<double, RandomNumberColumns> {
        return {dist.unif_rng(), dist.norm_rng(0, 1)};
    });
}

void BayesRBase::prepare(BayesRKernel *kernel)
{
    // Empty in BayesRBase
    (void) kernel; // Unused
}

void BayesRBase::readWithSharedLock(BayesRKernel *kernel)
{
    // Empty in BayesRBase
    (void) kernel; // Unused
}

void BayesRBase::writeWithUniqueLock(BayesRKernel *kernel)
{
    // Empty in BayesRBase
    (void) kernel; // Unused
}

int BayesRBase::runGibbs(AnalysisGraph *analysis)
{
    if (!analysis) {
        std::cout << "Cannot run Gibbs analysis without a flow graph!" << std::endl;
        return 1;
    }

    setAsynchronous(analysis->isAsynchronous());

    const unsigned int M(m_data->numSnps);
    const unsigned int N(m_data->numInds);
    const int K(int(m_cva.cols()) + 1);
    const unsigned int nGroups(m_data->numGroups);
    const unsigned int nF(m_opt->fixedEffectNumber);

    init(K, M, N);

    SampleWriter writer;
    writer.setFileName(m_outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.setFixedCount(nF);
    writer.openGroups(nGroups);

    LogWriter iterLogger;
    VectorXd  iterLog(10);

    if(m_showDebug)
    {

      iterLogger.setFileName(m_iterLogFile);
      iterLogger.open();
    }

    // Sampler variables
    VectorXd sample(3*M+3+nGroups+nF+N); // varible containg a sambple of all variables in the model, M marker effects, M component assigned to markers, M acum values, sigmaE, sigmaG, mu, iteration number and Explained variance
    std::vector<unsigned int> markerI(M);
    std::iota(markerI.begin(), markerI.end(), 0);

    //fixed effects vector & iterator
    m_gamma = VectorXd(nF);
    m_gamma.setZero();
    std::vector<unsigned int> xI(nF);
    std::iota(xI.begin(), xI.end(), 0);
    
    std::cout << "Number of groups: " << nGroups << std::endl
              << "Running Gibbs sampling" << endl;

    const auto t1 = std::chrono::high_resolution_clock::now();

    // This for MUST NOT BE PARALLELIZED, IT IS THE MARKOV CHAIN
    m_components = VectorXd::Zero(M);

    long meanIterationTime = 0;
    long meanFlowGraphIterationTime = 0;

    for (unsigned int iteration = 0; iteration < m_maxIterations; iteration++) {
        // Output progress
        const auto startTime = std::chrono::high_resolution_clock::now();
        //if (iteration > 0 && iteration % unsigned(std::ceil(max_iterations / 10)) == 0)
        std::cout << "iteration " << iteration << ": ";
        double old_mu=m_mu;

    // we delegate the Mu update to the descendents
        updateMu(old_mu,(double)N);
        const auto muTime = std::chrono::high_resolution_clock::now();
        prepareForAnylsis();

        std::random_shuffle(markerI.begin(), markerI.end());

        m_m0 = 0;
        m_v = MatrixXd::Zero(nGroups, K);
        m_betasqnG = VectorXd::Zero(nGroups);

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution.
        // The flow graph is constructed to allow the data to be decompressed in parallel for enforce sequential processing of each column
        // in turn. HOwever, within each column we make use of Intel TBB's parallel_for to parallelise the operations on the large vectors
        // of data.
        const auto flowGraphStartTime = std::chrono::high_resolution_clock::now();
        analysis->exec(this, N, M, markerI);
        const auto flowGraphEndTime = std::chrono::high_resolution_clock::now();
	
        // Fixed effects estimation
        // ---------------------
        double dNm1 = (double)(N - 1);
        if (nF>0) {
                std::random_shuffle(xI.begin(), xI.end());
                double gamma_old, num_f, denom_f;
                double sigE_sigF = m_sigmaE / m_sigmaF;

                for (int i=0; i<nF; i++) {
                    gamma_old = m_gamma(xI[i]);
                    num_f     = 0.0;
                    denom_f   = 0.0;                    
                    for (int k=0; k<N; k++){                  
                      num_f += m_data->X(k, xI[i]) * (m_epsilon[k] + gamma_old * m_data->X(k, xI[i]));
		    }
		      denom_f = dNm1 + sigE_sigF;
                      m_gamma(xI[i]) = m_dist.norm_rng(num_f/denom_f, m_sigmaE/denom_f);
                    
                    for (int k = 0; k<N ; k++) {
		      m_epsilon[k] = m_epsilon[k] + (gamma_old - m_gamma(xI[i])) * m_data->X(k, xI[i]);                        
                    }
                }
                m_sigmaF = s02F;
        }	

    const auto sEstartTime = std::chrono::high_resolution_clock::now();
        const double epsilonSqNorm = m_epsilon.squaredNorm();
        m_sigmaE = m_dist.inv_scaled_chisq_rng(m_v0E + N, (epsilonSqNorm + m_v0E * m_s02E) / (m_v0E + N));
    const auto sEendTime = std::chrono::high_resolution_clock::now();

        const auto sGstartTime = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < nGroups; i++) {
            m_m0 = m_v.row(i).sum() - m_v.row(i)(0);
            m_sigmaG[i] = m_dist.inv_scaled_chisq_rng(m_v0G + m_m0, (m_betasqnG(i) * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0));
            m_pi.row(i) = m_dist.dirichlet_rng(m_v.row(i).array() + 1.0);
        }
        const auto sGendTime = std::chrono::high_resolution_clock::now();

    if (iteration >= m_burnIn && iteration % m_thinning == 0) {
            sample << iteration, m_mu, m_beta, m_sigmaE, m_sigmaG, m_gamma, m_components, m_acum, m_epsilon;
            writer.write(sample);
        }

        const auto endTime = std::chrono::high_resolution_clock::now();
        const auto iterationDuration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
        const auto flowGraphDuration = std::chrono::duration_cast<std::chrono::milliseconds>(flowGraphEndTime - flowGraphStartTime).count();
        std::cout << static_cast<double>(iterationDuration) / 1000.0 << "s (" << static_cast<double>(flowGraphDuration) / 1000.0 << "s)" << std::endl;
        meanIterationTime += iterationDuration;
        meanFlowGraphIterationTime += flowGraphDuration;
    if (m_showDebug)
    {
      const auto muDuration = std::chrono::duration_cast<std::chrono::microseconds>(muTime - startTime).count();
      const auto sigmaGDuration = std::chrono::duration_cast<std::chrono::microseconds>(sGendTime - sGstartTime).count();
      const auto sigmaEDuration = std::chrono::duration_cast<std::chrono::microseconds>(sEendTime - sEstartTime).count();
      iterLog<<iteration,m_m0,m_sigmaG,m_sigmaE,m_mu
        ,static_cast<double>(muDuration)
        ,static_cast<double>(sigmaGDuration)
        ,static_cast<double>(sigmaEDuration)
        ,static_cast<double>(flowGraphDuration)
        ,static_cast<double>(iterationDuration);
      iterLogger.write(iterLog);
    }
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

void BayesRBase::processColumn(const KernelPtr &kernel)
{
    auto * bayesKernel = dynamic_cast<BayesRKernel*>(kernel.get());
    assert(bayesKernel);

    const unsigned int N(m_data->numInds);
    const double NM1 = double(N - 1);
    const int K(int(m_cva.cols()) + 1);
    const int km1 = K - 1;
    double acum = 0.0;
    double beta_old;
    const auto t1c = std::chrono::high_resolution_clock::now();

    beta_old = m_beta(bayesKernel->marker->i);
    const auto group = m_data->G[bayesKernel->marker->i];

    prepare(bayesKernel);
    readWithSharedLock(bayesKernel);

    // muk for the zeroth component=0
    m_muk[0] = 0.0;

    const double sigmaG = m_sigmaG[group];

    // set variable cVa
    m_cVa.segment(1, km1) = m_cva.row(group);
    m_cVaI.segment(1, km1) = m_cVa.segment(1, km1).cwiseInverse();

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / sigmaG;
    m_denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();

    const auto num_begin = std::chrono::high_resolution_clock::now();
    const double num = bayesKernel->computeNum(m_epsilon, beta_old);
    const auto num_end = std::chrono::high_resolution_clock::now();
    //The rest of the algorithm remains the same
     const auto beta_begin = std::chrono::high_resolution_clock::now();
    // muk for the other components is computed according to equaitons
    m_muk.segment(1, km1) = num / m_denom.array();

    // Update the log likelihood for each component
    VectorXd logL(K);
    const double logLScale = sigmaG / m_sigmaE * NM1;
    logL = m_pi.row(group).array().log(); // First component probabilities remain unchanged
    logL.segment(1, km1) = logL.segment(1, km1).array()
            - 0.5 * ((logLScale * m_cVa.segment(1, km1).array() + 1).array().log())
            + 0.5 * (m_muk.segment(1, km1).array() * num) / m_sigmaE;

    if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) {
        acum = 0;
    } else {
        acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
    }

    m_acum(bayesKernel->marker->i) = acum;
    
    const double p = m_randomNumbers.at(kernel->marker->i).at(PIndex);
    const double randomNorm = m_randomNumbers.at(kernel->marker->i).at(RandomNormIndex);

    for (int k = 0; k < K; k++) {
        if (p <= acum) {
            //if zeroth component
            if (k == 0) {
                m_beta(bayesKernel->marker->i) = 0;
            } else {
                m_beta(bayesKernel->marker->i) = randomNorm * std::sqrt(m_sigmaE/m_denom[k-1]) + m_muk[k];
                m_betasqnG(group) += pow(m_beta(bayesKernel->marker->i), 2);
            }
            m_v.row(group)(k)+=1.0;
            m_components[bayesKernel->marker->i] = k;
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
    const auto beta_end = std::chrono::high_resolution_clock::now();
    const double beta_new = m_beta(bayesKernel->marker->i);

    //until here
    //we skip update if old and new beta equals zero
    const bool skipUpdate = beta_old == 0.0 && beta_new == 0.0;

    const auto eps_begin = std::chrono::high_resolution_clock::now();
    if (!skipUpdate) {
        m_epsilon += *bayesKernel->calculateEpsilonChange(beta_old, beta_new);
        writeWithUniqueLock(bayesKernel);
    }
    const auto eps_end = std::chrono::high_resolution_clock::now();

    if(m_colLog)
    {
      Vector4d colLog;
      const auto numDuration = std::chrono::duration_cast<std::chrono::microseconds>(num_end - num_begin).count();
      const auto betaDuration = std::chrono::duration_cast<std::chrono::microseconds>(beta_end - beta_begin).count();
      const auto epsDuration = std::chrono::duration_cast<std::chrono::microseconds>(eps_end - eps_begin).count();
      colLog<< bayesKernel->marker->i, static_cast<double>(numDuration), static_cast<double>(betaDuration),static_cast<double>(epsDuration);
      m_colWriter.write(colLog);
    }
   // const auto durationc = std::chrono::duration_cast<std::chrono::microseconds>(t2c - t1c).count();

}

std::unique_ptr<AsyncResult> BayesRBase::processColumnAsync(const KernelPtr &kernel)
{
    auto * bayesKernel = dynamic_cast<BayesRKernel*>(kernel.get());
    assert(bayesKernel);

    const auto group = m_data->G[bayesKernel->marker->i];
    const double sigmaG = m_sigmaG[group];
    double component = 0;
    auto result = std::make_unique<AsyncResult>();
    // to keep track of the column processing time
    const auto t1c = std::chrono::high_resolution_clock::now();
    prepare(bayesKernel);

    // The elements of these members are only accessed by the thread we are in
    result->beta = m_beta(bayesKernel->marker->i);
    result->betaOld = result->beta;
    component = m_components(bayesKernel->marker->i);

    double num = 0;
    {
        // Use a shared lock to allow multiple threads to read updates
        std::shared_lock lock(m_mutex);
        readWithSharedLock(bayesKernel);//here we are reading the column and also epsilonsum
        num = bayesKernel->computeNum(m_epsilon, result->betaOld);
    }

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / sigmaG;

    const double NM1 = static_cast<double>(m_data->numInds) - 1.0;
    const int K(int(m_cva.cols()) + 1);
    const int km1 = K - 1;

    VectorXd cVa = VectorXd::Zero(K);
    cVa.segment(1, km1) = m_cva.row(group);
    VectorXd cVaI = VectorXd::Zero(K);
    cVaI.segment(1, km1) = cVa.segment(1, km1).cwiseInverse();

    VectorXd denom = NM1 + sigmaEOverSigmaG * cVaI.segment(1, km1).array();

    // muk for the zeroth component=0
    VectorXd muk(K);
    muk[0] = 0.0;
    // muk for the other components is computed according to equaitons
    muk.segment(1, km1) = num / denom.array();

    // Update the log likelihood for each component
    VectorXd logL(K);
    const double logLScale = sigmaG / m_sigmaE * NM1;
    logL = m_pi.row(group).array().log(); // First component probabilities remain unchanged
    logL.segment(1, km1) = logL.segment(1, km1).array()
            - 0.5 * ((logLScale * cVa.segment(1, km1).array() + 1).array().log())
            + 0.5 * (muk.segment(1, km1).array() * num) / m_sigmaE;


    double acum = 0.0;
    if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) {
        acum = 0;
    } else {
        acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
    }

    const double p = m_randomNumbers.at(kernel->marker->i).at(PIndex);
    const double randomNorm = m_randomNumbers.at(kernel->marker->i).at(RandomNormIndex);

    result->v = std::make_unique<VectorXd>(VectorXd::Zero(K));
    for (int k = 0; k < K; k++) {
        if (p <= acum) {
            //if zeroth component
            if (k == 0) {
                result->beta = 0;
            } else {
                result->beta = randomNorm * std::sqrt(m_sigmaE/denom[k-1]) + muk[k];
            }
            (*result->v)(k) += 1.0;
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

    // Only update m_epsilon if required
    const bool skipUpdate = result->betaOld == 0.0 && result->beta == 0.0;

    // Update our local copy of epsilon to minimise the amount of time we need to hold the unique lock for.
    if (!skipUpdate) {
          // this  also updates epsilonSum!
        result->deltaEpsilon = bayesKernel->calculateEpsilonChange(result->betaOld, result->beta);
        // now marker->epsilonSum now contains only delta_epsilonSum
    }

    // These updates do not need to be atomic
    m_beta(bayesKernel->marker->i) = result->beta;
    m_components(bayesKernel->marker->i) = component;

    //info on the running time of the column processing, would be very useful to have it as an option, and output it to a file
    //const auto t2c = std::chrono::high_resolution_clock::now();
   // const auto durationc = std::chrono::duration_cast<std::chrono::microseconds>(t2c - t1c).count();

   // cout<<"marker : "<< marker->i << " duration : " <<durationc<<endl;

    return result;
}

void BayesRBase::doThreadSafeUpdates(const ConstAsyncResultPtr &result)
{
    assert(result);

    // No mutex required here - thread_safe_update_node is serial, therefore
    // only one runs at any time. m_v is not accessed elsewhere whilst the
    // flow graph is running.
    m_v += *result->v;
}

void BayesRBase::updateGlobal(const KernelPtr& kernel,
                              const ConstAsyncResultPtr& result)
{
    assert(kernel);
    assert(result);

    std::unique_lock lock(m_mutex);
    m_epsilon += *result->deltaEpsilon;
    m_betasqnG[m_data->G[kernel->marker->i]] += pow(result->beta, 2);
}

void BayesRBase::printDebugInfo() const
{

}

/* I deferred this to its own class
void BayesRBase::printColumnDebugInfo() const
{
  cout<<"Marker, ": << endl;//TODO
  cout<<"marker_update"<<endl//TODO
  cout<<"num_update" <<endl;//TODO
  cout<<"deltaEps_update"<<<endl;//TODO
}

void BayesRBase::printGlobalDebugInfo()const
{
  cout << "Marker"<<endl;
  cout << "Epsilon_update"
  cout << "EpsilonSum_update"
  cout << "betasqn_update"
}
*/

