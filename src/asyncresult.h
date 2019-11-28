#ifndef ASYNCRESULT_H
#define ASYNCRESULT_H

#include <Eigen/Eigen>
#include <memory>

using namespace Eigen;

struct AsyncResult {
    double betaOld = 0.0;
    double beta = 0.0;
    std::unique_ptr<MatrixXd> v;
    double component = 0.0;

    std::unique_ptr<VectorXd> deltaEpsilon;

#if defined(EPSILON_TIMING_ENABLED)
    double count = 0.0;
#endif
};

#endif // ASYNCRESULT_H
