/*
 * BayesRRm.cpp
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#include "DenseBayesRRmz.hpp"
#include "limitsequencegraph.hpp"
#include "parallelgraph.hpp"

#include <mutex>

DenseBayesRRmz::DenseBayesRRmz(const Data *data, Options &opt)
    : BayesRBase (data, opt)
{
    if (opt.analysisType == "PPAsyncBayes") {
        m_flowGraph.reset(new ParallelGraph(this, opt.numThread));
        setAsynchronous(true);
    } else {
        m_flowGraph.reset(new LimitSequenceGraph(this, opt.numThread));
    }
}

DenseBayesRRmz::~DenseBayesRRmz()
{
}

void DenseBayesRRmz::processColumn(unsigned int marker, const Map<VectorXd> &Cx)
{
    const unsigned int N(m_data->numInds);
    const double NM1 = double(N - 1);
    const int K(int(m_cva.size()) + 1);
    const int km1 = K - 1;
    double acum = 0.0;
    double beta_old;

    beta_old = m_beta(marker);

    // Now y_tilde = Y-mu - X * beta + X.col(marker) * beta(marker)_old
    if (m_components(marker) != 0.0) {
        m_y_tilde = m_epsilon + beta_old * Cx;
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
      const double num = Cx.dot(m_y_tilde);
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
        m_epsilon = m_y_tilde - m_beta(marker) * Cx;
    } else {
        m_epsilon = m_y_tilde;
    }
    // Now epsilon contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new
}

std::tuple<double, double> DenseBayesRRmz::processColumnAsync(unsigned int marker, const Map<VectorXd> &Cx)
{
    // Lock and take local copies of needed variabls
    // [*] m_beta(marker) rwr - used, updated, then used - per column, could take a copy and update at end
    // [*] m_betasqn w - updated here, used in BayezRRmz::runGibbs
    // [*] m_components(marker) rwr - used, updated, then used - per column, could take a copy and update at end
    // [*] m_epsilon rw - used throughout, then updated, used in BayezRRmz::runGibbs
    // [*] m_v w - updated here, used in BayezRRmz::runGibbs

    // [*] m_dist r - the engine is not thread safe

    // Temporaries
    // - cost of locking vs allocating per iteration?
    // [*] m_denom wr - computed from m_cVaI
    // [*] m_muk wr - computed from m_cVaI

    // m_data.numInds r - could be a member?
    // m_cva.size() r - could be a member?
    // m_sigmaE r - calculated in BayesRRmz::init
    // m_sigmaG r - calculated in BayesRRmz::init, updated in BayezRRmz::runGibbs
    // m_pi r - calculated in BayesRRmz::init, updated in BayezRRmz::runGibbs
    // m_cVa r - calculated in BayesRRmz::init
    // m_cVaI r - calculated in BayesRRmz::init

    double beta = 0;
    double component = 0;
    VectorXd y_tilde(m_data->numInds);
    VectorXd epsilon(m_data->numInds);

    {
        // Use a shared lock to allow multiple threads to read updates
        std::shared_lock lock(m_mutex);

        // std::memcpy is faster than epsilon = m_epsilon which compiles down to a loop over pairs of
        // doubles and uses _mm_load_pd(source) SIMD intrinsics. Just be careful if we change the type
        // contained in the vector back to floats.
        std::memcpy(y_tilde.data(), m_async_epsilon.data(), static_cast<size_t>(epsilon.size()) * sizeof(double));
        beta = m_beta(marker);
        component = m_components(marker);
    }
    const double beta_old = beta;

    // Note that we assign y_tilde = m_epsilon above with the memcpy.
    // Now y_tilde = Y-mu - X * beta + X.col(marker) * beta(marker)_old
    if (component != 0.0)
        y_tilde += beta_old * Cx;

    // We compute the dot product to save computations
    const double num = Cx.dot(y_tilde);

    // Do work

    // We compute the denominator in the variance expression to save computations
    const double sigmaEOverSigmaG = m_sigmaE / m_sigmaG;

    const double NM1 = double(m_data->numInds - 1);
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
        y_tilde -= beta * Cx;
    }
    // Now y_tilde contains Y-mu - X*beta + X.col(marker) * beta(marker)_old - X.col(marker) * beta(marker)_new

    // Lock to write updates (at end, or perhaps as updates are computed)
    {
        // Use a unique lock to ensure only one thread can write updates
        std::unique_lock lock(m_mutex);
        if (!skipUpdate) {
            std::memcpy(m_async_epsilon.data(), y_tilde.data(), static_cast<size_t>(y_tilde.size()) * sizeof(double));
            m_betasqn += beta * beta - beta_old * beta_old;
        }
        m_v += v;
    }

    // These updates do not need to be atomic
    m_beta(marker) = beta;
    m_components(marker) = component;

    return {beta_old, beta};
}

void DenseBayesRRmz::updateGlobal(double beta_old, double beta, const Map<VectorXd> &Cx)
{
    // No mutex required here whilst m_globalComputeNode uses the serial policy
    m_epsilon -= Cx * (beta - beta_old);
}

void DenseBayesRRmz::init(int K, unsigned int markerCount, unsigned int individualCount)
{
    BayesRBase::init(K, markerCount, individualCount);

    m_async_epsilon = VectorXd(individualCount);
}

void DenseBayesRRmz::prepareForAnylsis()
{
    if (m_isAsync)
        std::memcpy(m_async_epsilon.data(), m_epsilon.data(), static_cast<size_t>(m_epsilon.size()) * sizeof(double));
}
