#ifndef BAYESWCOMMON_H_
#define BAYESWCOMMON_H_

#include <Eigen/Eigen>

using namespace Eigen;

// Two structures for ARS
struct pars{
        /* Common parameters for the densities */
        VectorXd epsilon;			// epsilon per subject (before each sampling, need to remove the effect of the sampled parameter and then carry on

        VectorXd mixture_classes; // Vector to store mixture component C_k values

        int used_mixture; //Write the index of the mixture we decide to use

        /* Store the current variables */
        double alpha, sigma_b;

        /* Beta_j - specific variables */
        VectorXd X_j;

        /*  of sum(X_j*failure) */
        double sum_failure;

        /* Mu-specific variables */
        double sigma_mu;
        /* sigma_b-specific variables */
        double alpha_sigma, beta_sigma;

        /* Number of events (sum of failure indicators) */
        double d;

        /* Help variable for storing sqrt(2sigma_b)	 */
        double sqrt_2sigmab;

};

struct pars_beta_sparse{
        /* Common parameters for the densities */
        VectorXd mixture_classes; // Vector to store mixture component C_k values

        int used_mixture; //Write the index of the mixture we decide to use

        /* Store the current variables */
        double alpha, sigma_b;

        /* Beta_j - specific variables */
        double vi_0, vi_1, vi_2; // Sums of vi elements

        // Mean, std dev and their ratio for snp j
        double mean, sd, mean_sd_ratio;

        /*  of sum(X_j*failure) */
        double sum_failure;
};

#endif /* BAYESWCOMMON_H_ */
