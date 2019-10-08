/*
 * BayesRRm.h
 *
 *  Created on: 5 Sep 2018
 *      Author: admin
 */

#ifndef SRC_BAYESRBASE_H_
#define SRC_BAYESRBASE_H_

#include "analysis.h"
#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"
#include "colwriter.h"

#include <Eigen/Eigen>
#include <memory>
#include <shared_mutex>

class AnalysisGraph;

struct BayesRKernel;

class BayesRBase : public Analysis
{
public:
    explicit BayesRBase(const Data *data, const Options *opt);

    int runGibbs(AnalysisGraph* analysis) override; // where we run Gibbs sampling over the parametrised model

    void processColumn(const KernelPtr &kernel) override;

    std::unique_ptr<AsyncResult> processColumnAsync(const KernelPtr &kernel) override;
    void doThreadSafeUpdates(const ConstAsyncResultPtr& result) override;
    void updateGlobal(const KernelPtr& kernel, const ConstAsyncResultPtr &result) override;

    virtual void updateMu(double old_mu, double N)=0;

    void setDebugEnabled(bool enabled) { m_showDebug = enabled; }
    bool isDebugEnabled() const { return m_showDebug; }

protected:
    const string        m_outputFile;
    const string        m_iterLogFile; //debug file for iteration quantities
    const unsigned int  m_seed;
    const unsigned int  m_maxIterations;
    const unsigned int  m_burnIn;
    const unsigned int  m_thinning;
    const double        m_sigma0=0.0001;
    const double        m_v0E;
    const double        m_s02E;
    const double        m_v0G;
    const double        m_s02G;
    MatrixXd            m_cva;
    Distributions_boost m_dist;
    bool m_showDebug;

    bool m_colLog=false;       //log for columns
    ColWriter m_colWriter;//writer for log for columns
    string m_colLogFile;

  // Component variables
    MatrixXd m_priorPi;   // prior probabilities for each component
    MatrixXd m_pi;        // mixture probabilities
    VectorXd m_cVa;       // component-specific variance
    VectorXd m_muk;       // mean of k-th component marker effect size
    VectorXd m_denom;     // temporal variable for computing the inflation of the effect variance for a given non-zero componnet
    int m_m0;             // total number of markers in model
    MatrixXd m_v;         // variable storing the component assignment
    VectorXd m_cVaI;      // inverse of the component variances

    // Mean and residual variables
    double m_mu;          // mean or intercept
    VectorXd m_sigmaG;    // genetic variance
    double m_sigmaE;      // residuals variance

    // Linear model variables
    VectorXd m_beta;       // effect sizes
    VectorXd m_acum;       // acum values, posterior inclusion probabilities per locus
    VectorXd m_y_tilde;    // variable containing the adjusted residuals to exclude the effects of a given marker
    VectorXd m_epsilon;    // variable containing the residuals
    VectorXd m_betasqnG;
    double m_epsilonSum=0.0;
    VectorXd m_y;
    VectorXd m_components;
    
    //fixed effects
    VectorXd m_gamma;     // fixed effects coefficients
    const double       s02F    = 1.0;
    double m_sigmaF;      // covariates variance if using ridge;

    static const std::size_t PIndex = 0;
    static const std::size_t RandomNormIndex = 1;
    constexpr static const std::size_t RandomNumberColumns = 2; // p, randomNorm
    std::vector<std::array<double, RandomNumberColumns>> m_randomNumbers;

    bool m_isAsync = false;

    mutable std::shared_mutex m_mutex;

    void setAsynchronous(bool async) { m_isAsync = async; }

    virtual void init(int K, unsigned int markerCount, unsigned int individualCount);

    virtual void prepareForAnylsis();

    virtual void prepare(BayesRKernel *kernel);
    virtual void readWithSharedLock(BayesRKernel *kernel);
    virtual void writeWithUniqueLock(BayesRKernel *kernel);

    void printDebugInfo() const;
};

#endif /* SRC_BAYESRBASE_H_ */
