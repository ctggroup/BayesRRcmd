#ifndef KERNEL_H
#define KERNEL_K

#include "marker.h"

#include <memory>

struct Kernel
{
    explicit Kernel(const Marker *marker);
    virtual ~Kernel();

    std::unique_ptr<const Marker> marker = nullptr;
};

#endif // KERNEL_H
