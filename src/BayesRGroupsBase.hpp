
#ifndef SRC_BAYESRGROUPSBASE_H_
#define SRC_BAYESRGROUPSBASE_H_

#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"
#include "BayesRBase.hpp"

#include <Eigen/Eigen>
#include <memory>
#include <shared_mutex>

class AnalysisGraph;

struct Marker;
class MarkerBuilder;

class BayesRGroupsBase : public BayesRBase
{
public:
    BayesRGroupsBase(const Data *m_data, Options &m_opt);
    virtual ~BayesRGroupsBase();

//  virtual MarkerBuilder* markerBuilder() const = 0;
//    virtual IndexEntry indexEntry(unsigned int i) const;
//    virtual bool compressed() const;
//    virtual unsigned char* compressedData() const;
//    virtual std::string preprocessedFile() const;

    int runGibbs(AnalysisGraph* analysis); // where we run Gibbs sampling over the parametrised model

    virtual void processColumn(Marker *marker) override;

//  virtual std::tuple<double, double,VectorXd> processColumnAsync(Marker *marker);
//  virtual void updateMu(double old_mu, double N)=0;
//  void setDebugEnabled(bool enabled) { m_showDebug = enabled; }
//  bool isDebugEnabled() const { return m_showDebug; }

protected:
//    const Data          *m_data; // data matrices
//    Options             &m_opt;
//    const string        m_bedFile; // bed file
//    const string        m_outputFile;
//    const unsigned int  m_seed;
//    const unsigned int  m_maxIterations;
//    const unsigned int  m_burnIn;
//    const unsigned int  m_thinning;
//    const double        m_sigma0  = 0.0001;
//    const double        m_v0E     = 0.0001;
//    const double        m_s02E    = 0.0001;
//    const double        m_v0G     = 0.0001;
//    const double        m_s02G    = 0.0001;

//    Distributions_boost m_dist;
//    bool m_usePreprocessedData;
//    bool m_showDebug;
    Eigen::MatrixXd     m_cva;

    // Component variables
    //TODO check if there is no problem by hiding the base class attributes
    MatrixXd m_priorPi;   // prior probabilities for each component
    MatrixXd m_pi;        // mixture probabilities
//    VectorXd m_cVa;       // component-specific variance
//    VectorXd m_muk;       // mean of k-th component marker effect size
//    VectorXd m_denom;     // temporal variable for computing the inflation of the effect variance for a given non-zero component
//    int m_m0;             // total number of markers in model
    MatrixXd m_v;         // variable storing the component assignment //TODO
//    VectorXd m_cVaI;      // inverse of the component variances

    // Mean and residual variables
//    double m_mu;          // mean or intercept
//    double m_sigmaG;      // genetic variance
//    double m_sigmaE;      // residuals variance
    VectorXd m_sigmaGG;	  // group specific genetic variance //TODO

    // Linear model variables
//    VectorXd m_beta;       // effect sizes
//    VectorXd m_y_tilde;    // variable containing the adjusted residuals to exclude the effects of a given marker
//    VectorXd m_epsilon;    // variable containing the residuals
//    double m_betasqn = 0.0;
//    double m_epsilonSum = 0.0;
//    VectorXd m_y;
//    VectorXd m_components;
    VectorXd m_betasqnG;	//TODO

//    bool m_isAsync = false;
//
//    VectorXd m_asyncEpsilon;
//
//    mutable std::shared_mutex m_mutex;
//    mutable std::mutex m_rngMutex;
//
//    void setAsynchronous(bool async) { m_isAsync = async; }

    virtual void init(int K, unsigned int markerCount, unsigned int individualCount, unsigned int groupCount);

//    virtual void prepareForAnylsis();
//
//    virtual void prepare(Marker *marker);
//    virtual void readWithSharedLock(Marker *marker);
//    virtual void writeWithUniqueLock(Marker *marker);
//
//    void printDebugInfo() const;
};

#endif /* SRC_BAYESRGROUPSBASE_H_ */
