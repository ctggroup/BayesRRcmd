#include "densemarker.h"

double DenseMarker::computeNum(VectorXd &epsilon, const double beta_old)
{
    if (component != 0.0)
        epsilon += beta_old * *Cx;

    return Cx->dot(epsilon);
}

void DenseMarker::updateEpsilon(VectorXd &epsilon, const double beta_old, const double beta)
{
    (void) beta_old; // Unused
    epsilon -= beta * *Cx;
}
