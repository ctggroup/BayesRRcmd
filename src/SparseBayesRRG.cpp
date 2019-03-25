#include "SparseBayesRRG.hpp"

#include "eigensparsedata.h"
#include "sparseparallelgraph.hpp"
#include "sparsesequentialanalysis.h"

#include <mutex>

// this code skeleton only highlights which would be the proposal of code changes
//Additional variables
// VectorXd means; //vector that contains the mean of each column of the bed file matrix
// VectorXd sds; //vector that contains the sd of each column of the bed file matrix
// VectorXd sqrdZ; //vector that contains the sum of squares of each column the bed file matrix
// VectorXd Zsum;  //vector that contains the sum of squares of each columnof the bed file matrix
// double epsilonSum //acumulator that updates the sum of epsilon vector


//SparseMatrix<double> Zg; //this is a sparse matrix that contains the uncentered and unscaled elements of the bed matrix
//OR
//std::vector<std:vector<int>> Zones(M)//vector containing the vectors the indexes of elements of the bed matrix which are one for each column
//std::vector<std:vector<int>> Ztwos(M)//vector containing the vectors the indexes of elements of the bed matrix which are two for each column
//

SparseBayesRRG::SparseBayesRRG(const EigenSparseData *data, Options &opt)
    : BayesRBase(data, opt)
    , m_sparseData(data)
    , m_means(data->means)
    , m_sds(data->sds)
    , m_sqrdZ(data->sqrdZ)
    , m_Zsum(data->Zsum)
{
    if (opt.analysisType == "PPAsyncBayes") {
        m_flowGraph.reset(new SparseParallelGraph(this, opt.numThread));
        setAsynchronous(true);
    } else {
        m_flowGraph.reset(new SparseSequentialAnalysis(this));
    }
}

SparseBayesRRG::~SparseBayesRRG()
{

}

void SparseBayesRRG::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    BayesRBase::init(K, markerCount, individualCount);

    m_asyncEpsilon = VectorXd(individualCount);

    m_epsilonSum = m_y.sum(); // we initialise with the current sum of y elements
    m_asyncEpsilonSum = m_epsilonSum;

    m_ones.setOnes(m_data->numInds);
}

void SparseBayesRRG::prepareForAnylsis()
{
    if (m_isAsync) {
        std::memcpy(m_asyncEpsilon.data(), m_epsilon.data(), static_cast<size_t>(m_epsilon.size()) * sizeof(double));
        m_asyncEpsilonSum = m_epsilonSum;
    }
}

void SparseBayesRRG::processColumn(unsigned int marker)
{
    const unsigned int N(m_data->numInds);
    const double NM1 = double(N - 1);
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    double acum = 0.0;
    double beta_old;

    beta_old = m_beta(marker);

    // muk for the zeroth component=0
    m_muk[0] = 0.0;

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;
    m_denom = NM1 + sigmaEOverSigmaG * m_cVaI.segment(1, km1).array();


    //OR the indexing solution, which using the development branch of eigen should be this
    //num = means(marker)*epsilonSum/sds(marker)+beta_old*sqrdZ(marker)-N*means(marker)/sds(marker) +(epsilon(Zones[marker]).sum()+2*epsilon(Ztwos[marker]).sum())/sds(marker)
    //maybe you can come up with a better way to index the elements of epsilon
    const double num = computeNum(marker, beta_old, m_epsilon,m_epsilonSum);

    //The rest of the algorithm remains the same
    
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
    const double beta_new = m_beta(marker);
    m_betasqn += beta_new * beta_new - beta_old * beta_old;
    
    //until here
    //we skip update if old and new beta equals zero
    const bool skipUpdate = beta_old == 0.0 && beta_new == 0.0;
    if (!skipUpdate) {
        const double dBeta = beta_old - beta_new;

        m_epsilon += computeEpsilonUpdate(marker, beta_old, beta_new);
        m_epsilonSum += computeEpsilonSumUpdate(marker, beta_old, beta_new);
    }
    // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
}

std::tuple<double, double> SparseBayesRRG::processColumnAsync(unsigned int marker)
{
    double beta = 0;
    double component = 0;
    VectorXd epsilon(m_data->numInds);
    double epsilonSum = 0;

    {
        // Use a shared lock to allow multiple threads to read updates
        std::shared_lock lock(m_mutex);

        // std::memcpy is faster than epsilon = m_epsilon which compiles down to a loop over pairs of
        // doubles and uses _mm_load_pd(source) SIMD intrinsics. Just be careful if we change the type
        // contained in the vector back to floats.
        std::memcpy(epsilon.data(), m_asyncEpsilon.data(), static_cast<size_t>(epsilon.size()) * sizeof(double));
        beta = m_beta(marker);
        component = m_components(marker);
        epsilonSum = m_asyncEpsilonSum;
    }
    const double beta_old = beta;
    const double N = static_cast<double>(m_data->numInds);
    const double NM1 = N - 1.0;

    const double num = computeNum(marker, beta_old, epsilon,epsilonSum);

    // Do work

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;

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
        epsilon += computeEpsilonUpdate(marker, beta_old, beta);
        epsilonSum += computeEpsilonSumUpdate(marker, beta_old, beta);
    }

    // Lock to write updates (at end, or perhaps as updates are computed)
    {
        // Use a unique lock to ensure only one thread can write updates
        std::unique_lock lock(m_mutex);
        if (!skipUpdate) {
            std::memcpy(m_asyncEpsilon.data(), epsilon.data(), static_cast<size_t>(epsilon.size()) * sizeof(double));
            m_asyncEpsilonSum = epsilonSum;
            m_betasqn += beta * beta - beta_old * beta_old;
        }
        m_v += v;
    }

    // These updates do not need to be atomic
    m_beta(marker) = beta;
    m_components(marker) = component;

    return {beta_old, beta};
}

void SparseBayesRRG::updateGlobal(const unsigned int marker, double beta_old, double beta)
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    m_epsilon += computeEpsilonUpdate(marker, beta_old, beta);
    m_epsilonSum += computeEpsilonSumUpdate(marker, beta_old, beta);
}

double SparseBayesRRG::computeNum(const unsigned int marker, const double beta_old, const VectorXd &epsilon,double epsilonSum) const
{
    //DANIEL here we either use the column of the sparse matrix or the two index vectors
    return beta_old * (static_cast<double>(m_data->numInds) - 1.0) - m_means(marker) * epsilonSum / m_sds(marker) + dot(marker, epsilon, m_sds(marker));;
    //OR the indexing solution, which using the development branch of eigen should be this
    //num = means(marker)*epsilonSum/sds(marker)+beta_old*sqrdZ(marker)-N*means(marker)/sds(marker) +(epsilon(Zones[marker]).sum()+2*epsilon(Ztwos[marker]).sum())/sds(marker)
}

double SparseBayesRRG::dot(const unsigned int marker, const VectorXd &epsilon, const double sd) const
{
    return m_sparseData->Zg[marker].dot(epsilon) / sd;
}

VectorXd SparseBayesRRG::computeEpsilonUpdate(const unsigned int marker, const double beta_old, const double beta) const
{
    const double dBeta = beta_old - beta;
    //Either
    return dBeta * m_sparseData->Zg[marker] / m_sds(marker) - dBeta * m_means(marker) / m_sds(marker) * m_ones;

    //OR
    //epsilon(Zones[marker])+=(beta_old-beta_new)/sds(marker)+(beta_new-beta_old)*means(marker)/sds(marker);
    //epsilon(Ztwos[marker])+=2*(beta_old-beta_new)/sds(marker)+(beta_new-beta_old)*means(marker)/sds(marker);
}

double SparseBayesRRG::computeEpsilonSumUpdate(const unsigned int marker, const double beta_old, const double beta) const
{
    //Regardless of which scheme, the update of epsilonSum is the same
    const double dBeta = beta_old - beta;
    return dBeta * m_Zsum(marker) / m_sds(marker) - dBeta * m_means(marker) * static_cast<double>(m_data->numInds) / m_sds(marker);
}
