/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "BayesRBase.hpp"
#include "samplewriter.h"
#include "analysisgraph.hpp"
#include "marker.h"
#include "logwriter.h"

#include <chrono>
#include <mutex>

BayesRBase::BayesRBase(const Data *data, Options &opt)
    : m_data(data)
    , m_opt(opt)
    , m_bedFile(opt.bedFile + ".bed")
    , m_outputFile(opt.mcmcSampleFile)
    , m_seed(opt.seed)
    , m_maxIterations(opt.chainLength)
    , m_burnIn(opt.burnin)
    , m_thinning(opt.thin)
    , m_dist(opt.seed)
    , m_usePreprocessedData(opt.analysisType == "PPBayes")
    , m_showDebug(opt.iterLog)
    , m_iterLogFile(opt.iterLogFile)
    , m_colLog(opt.colLog)
    , m_colLogFile(opt.colLogFile)
{
    assert(m_data);

    float* ptr =static_cast<float*>(&opt.S[0]);
    m_cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();
}

BayesRBase::~BayesRBase()
{
}

IndexEntry BayesRBase::indexEntry(unsigned int i) const
{
    if (!m_data)
        return {};

    return m_data->ppbedIndex[i];
}

bool BayesRBase::compressed() const
{
    return m_opt.compress;
}

unsigned char* BayesRBase::compressedData() const
{
    if (!m_data)
        return nullptr;

    return reinterpret_cast<unsigned char*>(m_data->ppBedMap);
}

std::string BayesRBase::preprocessedFile() const
{
    return ppFileForType(m_opt.dataType, m_opt.bedFile);
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
    m_epsilonSum=m_epsilon.sum();

    if(m_colLog)
    {
        m_colWriter.setFileName(m_colLogFile);
        m_colWriter.open();
    }
   
    
}

void BayesRBase::prepareForAnylsis()
{
    // Empty in BayesRBase
}

void BayesRBase::prepare(Marker *marker)
{
    // Empty in BayesRBase
    (void) marker; // Unused
}

void BayesRBase::readWithSharedLock(Marker *marker)
{
    // Empty in BayesRBase
    (void) marker; // Unused
}

void BayesRBase::writeWithUniqueLock(Marker *marker)
{
    // Empty in BayesRBase
    (void) marker; // Unused
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
    const int K(int(m_cva.size()) + 1);

    init(K, M, N);

    SampleWriter writer;
    writer.setFileName(m_outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.open();
    LogWriter iterLogger;
    VectorXd  iterLog(10);
    
    if(m_showDebug)
    {
      
      iterLogger.setFileName(m_iterLogFile);
      iterLogger.open();
    } 

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
        double old_mu=m_mu;

	// we delegate the Mu update to the descendents
        updateMu(old_mu,(double)N);
        const auto muTime = std::chrono::high_resolution_clock::now();
        prepareForAnylsis();

        std::random_shuffle(markerI.begin(), markerI.end());

        m_m0 = 0;
        m_v.setZero();

        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution.
        // The flow graph is constructed to allow the data to be decompressed in parallel for enforce sequential processing of each column
        // in turn. HOwever, within each column we make use of Intel TBB's parallel_for to parallelise the operations on the large vectors
        // of data.
        const auto flowGraphStartTime = std::chrono::high_resolution_clock::now();
        analysis->exec(this, N, M, markerI);
        const auto flowGraphEndTime = std::chrono::high_resolution_clock::now();

        m_m0 = int(M) - int(m_v[0]);
	const auto sGstartTime = std::chrono::high_resolution_clock::now();
        m_sigmaG = m_dist.inv_scaled_chisq_rng(m_v0G + m_m0, (m_betasqn * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0));
        const auto sGendTime = std::chrono::high_resolution_clock::now();

	const auto sEstartTime = std::chrono::high_resolution_clock::now();
        const double epsilonSqNorm = m_epsilon.squaredNorm();
        m_sigmaE = m_dist.inv_scaled_chisq_rng(m_v0E + N, (epsilonSqNorm + m_v0E * m_s02E) / (m_v0E + N));
	const auto sEendTime = std::chrono::high_resolution_clock::now();

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

void BayesRBase::processColumn(Marker *marker)
{
    const unsigned int N(m_data->numInds);
    const double NM1 = double(N - 1);
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    double acum = 0.0;
    double beta_old;
    const auto t1c = std::chrono::high_resolution_clock::now();

    beta_old = m_beta(marker->i);

    prepare(marker);
    readWithSharedLock(marker);

    // muk for the zeroth component=0
    m_muk[0] = 0.0;

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;
    m_denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();

    const auto num_begin = std::chrono::high_resolution_clock::now();
    const double num = marker->computeNum(m_epsilon, beta_old);
    const auto num_end = std::chrono::high_resolution_clock::now();
    //The rest of the algorithm remains the same
     const auto beta_begin = std::chrono::high_resolution_clock::now();
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
                m_beta(marker->i) = 0;
            } else {
                m_beta(marker->i) = m_dist.norm_rng(m_muk[k], m_sigmaE/m_denom[k-1]);
            }
            m_v[k] += 1.0;
            m_components[marker->i] = k;
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
    const double beta_new = m_beta(marker->i);
    m_betasqn += beta_new * beta_new - beta_old * beta_old;

    //until here
    //we skip update if old and new beta equals zero
    const bool skipUpdate = beta_old == 0.0 && beta_new == 0.0;

    const auto eps_begin = std::chrono::high_resolution_clock::now();
    if (!skipUpdate) {
      
        marker->updateEpsilon(m_epsilon, beta_old, beta_new);
        writeWithUniqueLock(marker);

    }
    const auto eps_end = std::chrono::high_resolution_clock::now();

    if(m_colLog)
    {
      Vector4d colLog;
      const auto numDuration = std::chrono::duration_cast<std::chrono::microseconds>(num_end - num_begin).count();
      const auto betaDuration = std::chrono::duration_cast<std::chrono::microseconds>(beta_end - beta_begin).count();
      const auto epsDuration = std::chrono::duration_cast<std::chrono::microseconds>(eps_end - eps_begin).count();
      colLog<< marker->i, static_cast<double>(numDuration), static_cast<double>(betaDuration),static_cast<double>(epsDuration);
      m_colWriter.write(colLog);
    }
   // const auto durationc = std::chrono::duration_cast<std::chrono::microseconds>(t2c - t1c).count();
   
}

std::tuple<double, double,VectorXd> BayesRBase::processColumnAsync(Marker *marker)
{
    double beta = 0;
    double component = 0;
    VectorXd epsilon(m_data->numInds);
    VectorXd deltaEps(m_data->numInds); //vector that will contain the delta epsilon message
    deltaEps.setZero(); // deltaEps=0
    // to keep track of the column processing time     
    const auto t1c = std::chrono::high_resolution_clock::now();
    prepare(marker);

    {
        // Use a shared lock to allow multiple threads to read updates
        std::shared_lock lock(m_mutex);

        // std::memcpy is faster than epsilon = m_epsilon which compiles down to a loop over pairs of
        // doubles and uses _mm_load_pd(source) SIMD intrinsics. Just be careful if we change the type
        // contained in the vector back to floats.
        // we copy global into local
        std::memcpy(epsilon.data(), m_epsilon.data(), static_cast<size_t>(epsilon.size()) * sizeof(double));
        beta = m_beta(marker->i);
        component = m_components(marker->i);
	 readWithSharedLock(marker);//here we are reading the column and also epsilonsum
    }
  
    const double beta_old = beta;

    const double num = marker->computeNum(epsilon, beta_old);

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;

    const double NM1 = static_cast<double>(m_data->numInds) - 1.0;
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    VectorXd denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();

    // muk for the zeroth component=0
    VectorXd muk(K);
    muk[0] = 0.0;
    // muk for the other components is computed according to equaitons
    muk.segment(1, km1) = num / denom.array();

    // Update the log likelihood for each component
    VectorXd logL(K);
    const double logLScale = m_sigmaG / m_sigmaE * NM1;
    logL = m_pi.array().log(); // First component probabilities remain unchanged
    logL.segment(1, km1) = logL.segment(1, km1).array()
            - 0.5 * ((logLScale * m_cVa.segment(1, km1).array() + 1).array().log())
            + 0.5 * (muk.segment(1, km1).array() * num) / m_sigmaE;


    double acum = 0.0;
    if (((logL.segment(1, km1).array() - logL[0]).abs().array() > 700).any()) {
        acum = 0;
    } else {
        acum = 1.0 / ((logL.array() - logL[0]).exp().sum());
    }

    double p = 0;
    std::vector<double> randomNumbers(static_cast<std::vector<double>::size_type>(K), 0);
    {
        // Generate all the numbers we are going to need in one go.
        // Use a unique lock to ensure only one thread can use the random number engine
        // at a time.
        std::unique_lock lock(m_rngMutex);
        p = m_dist.unif_rng();

        auto beginItr = randomNumbers.begin();
        std::advance(beginItr, 1);
        std::generate(beginItr, randomNumbers.end(), [&, k = 0] () mutable {
            ++k;
            return m_dist.norm_rng(muk[k], m_sigmaE/denom[k-1]);
        });
    }
    VectorXd v = VectorXd(K);
    v.setZero();
    for (int k = 0; k < K; k++) {
        if (p <= acum) {
            //if zeroth component
            if (k == 0) {
                beta = 0;
            } else {
                beta = randomNumbers.at(static_cast<std::vector<double>::size_type>(k));
            }
            v[k] += 1.0;
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
    const bool skipUpdate = beta_old == 0.0 && beta == 0.0;

    // Update our local copy of epsilon to minimise the amount of time we need to hold the unique lock for.
    if (!skipUpdate) {    
          // this  also updates epsilonSum!
        marker->updateEpsilon(deltaEps, beta_old, beta);    
        // now marker->epsilonSum now contains only delta_epsilonSum    
    }
    // In the new version of Async we do not synchronise epsilon Async, we will handle this through the global node
    // Lock to write updates (at end, or perhaps as updates are computed)
    
     {
      std::unique_lock lock(m_mutex);
      // Use a unique lock to ensure only one thread can write updates
         m_v += v; //maybe we can move this to the message struct
      }

    // These updates do not need to be atomic
    m_beta(marker->i) = beta;
    m_components(marker->i) = component;
    
    //info on the running time of the column processing, would be very useful to have it as an option, and output it to a file
    //const auto t2c = std::chrono::high_resolution_clock::now();
   // const auto durationc = std::chrono::duration_cast<std::chrono::microseconds>(t2c - t1c).count();
   
   // cout<<"marker : "<< marker->i << " duration : " <<durationc<<endl;

    return {beta_old, beta, deltaEps};
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
