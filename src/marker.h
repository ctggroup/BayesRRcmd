#ifndef MARKER_H
#define MARKER_H

#include <Eigen/Eigen>
#include <memory>

using namespace Eigen;

struct Marker
{
    virtual ~Marker();

    unsigned int i = 0;
    unsigned int numInds = 0;

    virtual double computeNum(VectorXd &epsilon,
                              const double beta_old) = 0;

    virtual void updateEpsilon(VectorXd &epsilon,
                               const double beta_old,
                               const double beta) = 0;
};

struct CompressedMarker
{
    std::shared_ptr<unsigned char[]> buffer = nullptr;
    unsigned long size = 0;
};

template<typename MarkerType>
CompressedMarker compress(const MarkerType* marker);

template<typename MarkerType>
void write(const MarkerType* marker, std::ostream *outStream);

#endif // MARKER_H
