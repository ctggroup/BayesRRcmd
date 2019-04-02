#ifndef MARKER_H
#define MARKER_H

#include <Eigen/Eigen>

using namespace Eigen;

struct Marker
{
    virtual ~Marker();

    unsigned int i = 0;

    virtual double computeNum(VectorXd &epsilon,
                              const double beta_old) = 0;

    virtual void updateEpsilon(VectorXd &epsilon,
                               const double beta_old,
                               const double beta) = 0;
};

template<typename MarkerType>
unsigned long compress(const MarkerType* marker,
                       unsigned char *outputBuffer,
                       unsigned long outputSize);

template<typename MarkerType>
void write(const MarkerType* marker, std::ofstream *outStream);

#endif // MARKER_H
