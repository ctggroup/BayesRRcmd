
#include "BayesRGroupsBase.hpp"
#include "samplewriter.h"
#include "analysisgraph.hpp"
#include "marker.h"

#include <chrono>
#include <mutex>

BayesRGroupsBase::BayesRGroupsBase(const Data *data, Options &opt): BayesRBase(data, opt)
{
    assert(m_data);

    //float* ptr =static_cast<float*>(&opt.S[0]);
    //m_cva = (Eigen::Map<Eigen::VectorXf>(ptr, static_cast<long>(opt.S.size()))).cast<double>();

    m_cva = MatrixXd(m_data->mS);
}

BayesRGroupsBase::~BayesRGroupsBase()
{
}

void BayesRGroupsBase::init(int K, unsigned int markerCount, unsigned int individualCount, unsigned int groupCount)
{
    // Component variables
	m_priorPi = MatrixXd(groupCount,K);		// prior probabilities for each component
    m_pi = MatrixXd(groupCount,K);          // mixture probabilities
    m_cVa = VectorXd(K);          			// component-specific variance
    m_muk = VectorXd (K);         			// mean of k-th component marker effect size
    m_denom = VectorXd(K - 1);    			// temporal variable for computing the inflation of the effect variance for a given non-zero componnet
    m_m0 = 0;                     			// total num ber of markes in model
    m_v = MatrixXd(groupCount,K);           // variable storing the component assignment
    m_cVaI = VectorXd(K);         			// inverse of the component variances

    // Mean and residual variables
    m_mu = 0.0;       						// mean or intercept
    m_sigmaG = 0.0;   						// genetic variance
    m_sigmaE = 0.0;   						// residuals variance
    m_sigmaGG = VectorXd(groupCount);

    // Linear model variables
    m_beta = VectorXd(markerCount);           // effect sizes
    m_y_tilde = VectorXd(individualCount);    // variable containing the adjusted residuals to exclude the effects of a given marker
    m_epsilon = VectorXd(individualCount);    // variable containing the residuals

    m_y = VectorXd();
    //Cx = VectorXd();

    m_betasqnG = VectorXd(groupCount);

    // Init the working variables
    const int km1 = K - 1;
    m_cVa[0] = 0;
    m_cVaI[0] = 0;

    for(int i=0; i < groupCount; i++){
        	m_priorPi.row(i)(0)=0.5;
        	for(int k=1;k<K;k++){
        		m_priorPi.row(i)(k)=0.5/K;
        	}
        }

    m_y_tilde.setZero();
    m_beta.setZero();

   	for(int i=0; i<groupCount;i++)
   		m_sigmaGG[i] = m_dist.beta_rng(1,1);

    m_pi = m_priorPi;

    m_y = (m_data->y.cast<double>().array() - m_data->y.cast<double>().mean());
    m_y /= sqrt(m_y.squaredNorm() / (double(individualCount - 1)));

    m_epsilon = (m_y).array() - m_mu;
    m_sigmaE = m_epsilon.squaredNorm() / individualCount * 0.5;
    m_epsilonSum=m_epsilon.sum();
}

int BayesRGroupsBase::runGibbs(AnalysisGraph *analysis)
{
    if (!analysis) {
        std::cout << "Cannot run Gibbs analysis without a flow graph!" << std::endl;
        return 1;
    }

    setAsynchronous(analysis->isAsynchronous());

    const unsigned int M(m_data->numSnps);
    const unsigned int N(m_data->numInds);
    const unsigned int nGroups(m_data->numGroups);
    const int K(m_data->mS.cols() + 1);

    init(K, M, N, nGroups);

    SampleWriter writer;
    writer.setFileName(m_outputFile);
    writer.setMarkerCount(M);
    writer.setIndividualCount(N);
    writer.openGroups(nGroups);

    // Sampler variables
    VectorXd sample(2*M+3+nGroups+N); // varible containg a sample of all variables in the model, M marker effects, M component assigned to markers, sigmaE, sigmaG, mu, iteration number and Explained variance
    std::vector<unsigned int> markerI(M);
    std::iota(markerI.begin(), markerI.end(), 0);

    std::cout << "Number of groups: " << nGroups << std::endl;

    std::cout << "Running Gibbs sampling" << std::endl;

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

        prepareForAnylsis();

        std::random_shuffle(markerI.begin(), markerI.end());

        m_m0 = 0;
        m_v.setZero();
        m_betasqnG.setZero();


        // This for should not be parallelized, resulting chain would not be ergodic, still, some times it may converge to the correct solution.
        // The flow graph is constructed to allow the data to be decompressed in parallel for enforce sequential processing of each column
        // in turn. HOwever, within each column we make use of Intel TBB's parallel_for to parallelise the operations on the large vectors
        // of data.

        const auto flowGraphStartTime = std::chrono::high_resolution_clock::now();
        analysis->exec(this, N, M, markerI);
        const auto flowGraphEndTime = std::chrono::high_resolution_clock::now();

        if (m_showDebug)
            printDebugInfo();
        const double epsilonSqNorm = m_epsilon.squaredNorm();
        m_sigmaE = m_dist.inv_scaled_chisq_rng(m_v0E + N, (epsilonSqNorm + m_v0E * m_s02E) / (m_v0E + N));

    	for(int i = 0; i < nGroups; i++){
    		m_m0 = m_v.row(i).sum() - m_v.row(i)(0);
    		m_sigmaGG[i] = m_dist.inv_scaled_chisq_rng(m_v0G + m_m0, (m_betasqnG(i) * m_m0 + m_v0G * m_s02G) / (m_v0G + m_m0));
    		m_pi.row(i) = m_dist.dirichilet_rng(m_v.row(i).array() + 1.0);
    	}

    	m_betasqn= m_betasqnG.sum();


        if (iteration >= m_burnIn && iteration % m_thinning == 0) {
            sample << iteration, m_mu, m_beta, m_sigmaE, m_sigmaGG, m_components, m_epsilon;
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

void BayesRGroupsBase::processColumn(Marker *marker)
{
    const unsigned int N(m_data->numInds);
    const double NM1 = double(N - 1);
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    double acum = 0.0;
    double beta_old;
    const auto t1c = std::chrono::high_resolution_clock::now();
    double sigmaG_process;
    beta_old = m_beta(marker->i);

    prepare(marker);
    readWithSharedLock(marker);

    // muk for the zeroth component=0
    m_muk[0] = 0.0;

    // get sigmaGG
    sigmaG_process = m_sigmaGG[m_data->G(marker->i)];

    // set variable cVa
    m_cVa.segment(1, km1) = m_cva.row(m_data->G(marker->i));
    m_cVaI.segment(1, km1) = m_cVa.segment(1, km1).cwiseInverse();

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / sigmaG_process;
    m_denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();

    const double num = marker->computeNum(m_epsilon, beta_old);

    //The rest of the algorithm remains the same

    // muk for the other components is computed according to equaitons
    m_muk.segment(1, km1) = num / m_denom.array();

    // Update the log likelihood for each component
    VectorXd logL(K);
    const double logLScale = sigmaG_process / m_sigmaE * NM1;
    logL = m_pi.row(m_data->G(marker->i)).array().log();
    // First component probabilities remain unchanged
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
                m_betasqnG(m_data->G(marker->i))+= pow(m_beta(marker->i),2);
            }
            m_v.row(m_data->G(marker->i))(k)+=1.0;
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
    const double beta_new = m_beta(marker->i);

    //until here
    //we skip update if old and new beta equals zero
    const bool skipUpdate = beta_old == 0.0 && beta_new == 0.0;
    if (!skipUpdate) {
        marker->updateEpsilon(m_epsilon, beta_old, beta_new);
        writeWithUniqueLock(marker);
    }
   //const auto t2c = std::chrono::high_resolution_clock::now();
   // const auto durationc = std::chrono::duration_cast<std::chrono::microseconds>(t2c - t1c).count();
   // cout<<"marker : "<< marker->i << " duration : " <<durationc<<endl;
}

