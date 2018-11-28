/*
 * BayesW.hpp
 *
 *  Created on: 26 Nov 2018
 *      Author: admin
 */

#ifndef BAYESW_HPP_
#define BAYESW_HPP_

#include "data.hpp"
#include "options.hpp"
#include "distributions_boost.hpp"

#include <Eigen/Eigen>

class BayesW
{
    Data            &data; // data matrices
    Options         &opt;
    const string    bedFile; // bed file
    const long      memPageSize; // size of memory
    const string    outputFile;
    const int       seed;
    const int       max_iterations;
    const int		burn_in;
    const int       thinning;
    const double	alpha_0  = 0.01;
    const double	kappa_0     = 0.01;
    const double    sigma_mu    = 100;
    const double    alpha_sigma  = 0.01;
    const double    beta_sigma   = 0.01;
  //  Eigen::VectorXd cva;
    Distributions_boost dist;

public:
    BayesW(Data &data, Options &opt, const long memPageSize);
    virtual ~BayesW();
    int runGibbs(); // where we run Gibbs sampling over the parametrised model
};



#endif /* BAYESW_HPP_ */
