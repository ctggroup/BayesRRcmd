#ifndef PARALLELALGO_H
#define PARALLELALGO_H

#include <Eigen/Eigen>

using namespace Eigen;

class Distributions_boost;

// Outer loop helpers
double parallelStepAndSumEpsilon(VectorXd &epsilon, double mu);
void parallelStepMuEpsilon(double &mu, VectorXd &epsilon, double sigmaEpsilon, double N, double sigmaE, Distributions_boost &dist);
double parallelSquaredNorm(const VectorXd &epsilon);

// Inner loop helpers
void parallelCastDouble(const VectorXf &in, VectorXd &out);
void parallelUpdateYTilde(VectorXd &y_tilde, const VectorXd &epsilon, const VectorXf &Cx, double beta);
void parallelUpdateEpsilon(VectorXd &epsilon, const VectorXd &y_tilde, const VectorXf &Cx, double beta);
double parallelDotProduct(const VectorXf &Cx, const VectorXd &y_tilde);
double parallelcwiseProduct(const VectorXf &Cx, const VectorXd &y_tilde);

#endif // PRALLELALGO_H
