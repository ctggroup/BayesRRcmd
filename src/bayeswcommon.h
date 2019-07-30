#ifndef BAYESWCOMMON_H_
#define BAYESWCOMMON_H_

#include <Eigen/Eigen>

using namespace Eigen;

struct pars_beta_sparse{
        /* Beta_j - specific variables */
        double vi_0, vi_1, vi_2; // Sums of vi elements

        // Mean, std dev and their ratio for snp j
        double mean, sd, mean_sd_ratio;

        /*  of sum(X_j*failure) */
        double sum_failure;
};

#endif /* BAYESWCOMMON_H_ */
