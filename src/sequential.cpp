#include "sequential.h"

#include "analysis.h"
#include "kernel.h"
#include "markercache.h"

#include <algorithm>
#include <iostream>

void Sequential::exec(Analysis *analysis,
                      unsigned int numInds,
                      unsigned int numSnps,
                      const std::vector<unsigned int> &markerIndices)
{
    (void) numInds; // Unused
    (void) numSnps; // Unused

    if (!analysis) {
        std::cerr << "Cannot run Sequential without bayes" << std::endl;
        return;
    }

    std::for_each(markerIndices.cbegin(), markerIndices.cend(), [&analysis](unsigned int i) {
        KernelPtr kernel = analysis->kernelForMarker(markerCache()->marker(i));
        analysis->processColumn(kernel);
    });
}
