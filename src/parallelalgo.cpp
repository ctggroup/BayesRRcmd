#include "parallelalgo.h"
#include "distributions_boost.hpp"
#include <tbb/tbb.h>

const size_t grainSize = 2000;


using namespace tbb;

double parallelStepAndSumEpsilon(VectorXd &epsilon, double mu)
{
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(epsilon.size());

    auto apply = [&](const blocked_range<size_t>& r, double initialValue) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());

        // Update epsilon. We substract previous value
        epsilon.segment(start, count) = epsilon.segment(start, count).array() + mu;
        const auto sum = initialValue + epsilon.segment(start, count).sum();
        return sum;
    };

    auto combine = [](double a, double b) { return a + b; };

    return parallel_reduce(blocked_range<unsigned long>(startIndex, endIndex, grainSize), 0.0,
                           apply, combine);
}

void parallelStepMuEpsilon(double &mu, VectorXd &epsilon, double sigmaEpsilon, double N, double sigmaE, Distributions_boost &dist)
{
    const double sigmaEpsilonOverN = sigmaEpsilon / N;
    const double sigmaEOverN = sigmaE / N;
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(epsilon.size());
    // Update mu
     mu = dist.norm_rng(sigmaEpsilonOverN, sigmaEOverN);

    parallel_for(blocked_range<unsigned long>(startIndex, endIndex, grainSize), [&](const blocked_range<size_t>& r) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());
        // We substract again now epsilon =Y-mu-X*beta
        epsilon.segment(start, count) = epsilon.segment(start, count).array() - mu;
    } );
}

double parallelSquaredNorm(const VectorXd &epsilon)
{
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(epsilon.size());

    auto apply = [&](const blocked_range<size_t>& r, double initialValue) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());

        const auto sum = initialValue + epsilon.segment(start, count).squaredNorm();
        return sum;
    };

    auto combine = [](double a, double b) { return a + b; };

    return parallel_reduce(blocked_range<unsigned long>(startIndex, endIndex, grainSize), 0.0,
                           apply, combine);
}

void parallelCastDouble(const VectorXf &in, VectorXd &out)
{
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(in.size());

    parallel_for(blocked_range<unsigned long>(startIndex, endIndex, grainSize), [&](const blocked_range<size_t>& r) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());
        out.segment(start, count) = in.segment(start, count).cast<double>();
    } );
}

void parallelUpdateYTilde(VectorXd &y_tilde, const VectorXd &epsilon, const VectorXf &Cx, double beta)
{
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(y_tilde.size());

    parallel_for(blocked_range<unsigned long>(startIndex, endIndex, grainSize), [&](const blocked_range<size_t>& r) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());
        y_tilde.segment(start, count) = epsilon.segment(start, count) + Cx.segment(start, count).cast<double>() * beta;
    } );
}

void parallelUpdateEpsilon(VectorXd &epsilon, const VectorXd &y_tilde, const VectorXf &Cx, double beta)
{
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(y_tilde.size());

    parallel_for(blocked_range<unsigned long>(startIndex, endIndex, grainSize), [&](const blocked_range<size_t>& r) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());
        epsilon.segment(start, count) = y_tilde.segment(start, count) - Cx.segment(start, count).cast<double>() * beta;
    } );
}

double parallelDotProduct(const VectorXf &Cx, const VectorXd &y_tilde)
{
    const unsigned long startIndex = 0;
    const unsigned long endIndex = static_cast<unsigned long>(y_tilde.size());

    auto apply = [&](const blocked_range<size_t>& r, double initialValue) {
        const long start = static_cast<long>(r.begin());
        const long count = static_cast<long>(r.end() - r.begin());
        const auto sum = initialValue + Cx.segment(start, count).cast<double>().dot(y_tilde.segment(start, count));
        return sum;
    };

    auto combine = [](double a, double b) { return a + b; };

    return parallel_reduce(blocked_range<unsigned long>(startIndex, endIndex, grainSize), 0.0,
                           apply, combine);
}
