#ifndef KERNEL_H
#define KERNEL_K

#include "marker.h"

#include <memory>

struct Kernel
{
    explicit Kernel(const ConstMarkerPtr &marker);
    virtual ~Kernel();

    ConstMarkerPtr marker = nullptr;

    virtual VectorXdPtr calculateEpsilonChange(const double beta_old,
                                               const double beta) = 0;
};

#endif // KERNEL_H
