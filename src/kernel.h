#ifndef KERNEL_H
#define KERNEL_K

#include "marker.h"

#include <memory>

struct Kernel
{
    explicit Kernel(const ConstMarkerPtr &marker);
    virtual ~Kernel();

    ConstMarkerPtr marker = nullptr;
};

#endif // KERNEL_H
