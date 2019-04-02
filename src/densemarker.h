#ifndef DENSEMARKER_H
#define DENSEMARKER_H

#include "marker.h"

struct DenseMarker : public Marker
{
    double component = 0;
    const Map<VectorXd> *Cx = nullptr;

    double computeNum(VectorXd &epsilon, const double beta_old) override;
    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;
};

#endif // DENSEMARKER_H
