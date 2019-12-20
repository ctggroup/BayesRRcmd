#ifndef ASYNCRESULT_H
#define ASYNCRESULT_H

#include <Eigen/Eigen>
#include <memory>

using namespace Eigen;

struct AsyncResult {
    // copied between hosts
    double betaOld = 0.0;
    double beta = 0.0;
    double component = 0.0;

    // local results
    std::unique_ptr<VectorXd> deltaEpsilon;
    std::unique_ptr<MatrixXd> v;

#if defined(EPSILON_TIMING_ENABLED)
    double count = 0.0;
#endif
};

#endif // ASYNCRESULT_H
