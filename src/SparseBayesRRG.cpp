#include "SparseBayesRRG.hpp"

#include "sparsedata.h"
#include "sparsesequentialanalysis.h"

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

SparseBayesRRG::SparseBayesRRG(const SparseData *data, Options &opt)
    : BayesRBase(data, opt)
    , m_sparseData(data)
    , m_means(data->means)
    , m_sds(data->sds)
    , m_sqrdZ(data->sqrdZ)
    , m_Zsum(data->Zsum)
{
    m_flowGraph.reset(new SparseSequentialAnalysis(this));
}

SparseBayesRRG::~SparseBayesRRG()
{

}

void SparseBayesRRG::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    BayesRBase::init(K, markerCount, individualCount);

    m_epsilonSum = m_y.sum(); // we initialise with the current sum of y elements
    m_ones.setOnes(m_data->numInds);
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

    const double num = computeNum(marker, beta_old, m_epsilon);

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
        m_epsilon += computeEpsilonUpdate(marker, beta_old, beta_new);
        m_epsilonSum += computeEpsilonSumUpdate(marker, beta_old, beta_new);
    }
    // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
}

double SparseBayesRRG::computeNum(const unsigned int marker, const double beta_old, const VectorXd &epsilon) const
{
    //DANIEL here we either use the column of the sparse matrix or the two index vectors
    return m_means(marker) * m_epsilonSum / m_sds(marker) + beta_old * (static_cast<double>(m_data->numInds) - 1.0) + dot(marker, epsilon, m_sds(marker));
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

#if 0
void SparseBayesRRG::updateGlobal(double beta_old, double beta, const Map<VectorXd> &Cx)
{
    // Not yet implemented
    (void)beta_old;
    (void)beta;
    (void)Cx;
    assert(false);

    // No mutex required here whilst m_globalComputeNode uses the serial policy
    //Either
    
    
    //OR
    //epsilon(Zones[marker])+=(beta_old-beta_new)/sds(marker)+(beta_new-beta_old)*means(marker)/sds(marker);
    //epsilon(Ztwos[marker])+=2*(beta_old-beta_new)/sds(marker)+(beta_new-beta_old)*means(marker)/sds(marker);

    //Regardless of which scheme, the update of epsilonSum is the same
    //    epsilonSum+= (beta_old-beta_new)*Zsum(marker)
    
}
#endif
