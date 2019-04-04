#ifndef MARKER_H
#define MARKER_H

#include <Eigen/Eigen>
#include <memory>

using namespace Eigen;

struct IndexEntry;

struct CompressedMarker
{
    std::shared_ptr<unsigned char[]> buffer = nullptr;
    unsigned long size = 0;
};

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

    virtual CompressedMarker compress() const = 0;
    virtual void decompress(unsigned char* data,
                            const IndexEntry& index) = 0;

    virtual size_t size() const { return 0; }
    virtual void write(std::ostream *outStream) const = 0;
};

#endif // MARKER_H
