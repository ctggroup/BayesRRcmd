#ifndef SPARSEMARKER_H
#define SPARSEMARKER_H

#include "marker.h"

struct SparseMarker : public Marker
{
    double mean = 0;
    double sd = 0;
    double sqrdZ= 0;
    double Zsum = 0;

    double epsilonSum = 0;

    double computeNum(VectorXd &epsilon, const double beta_old) override;

    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;

protected:
    virtual double computeNum(VectorXd &epsilon,
                              const double beta_old,
                              const double epsilonSum);

    virtual double dot(const VectorXd &epsilon) const = 0;

    virtual double computeEpsilonSumUpdate(const double beta_old,
                                           const double beta) const;
};

void updateStatistics(SparseMarker* marker,
                      unsigned int allele1,
                      unsigned int allele2);

std::size_t statisticsSize();

bool writeStatistics(const SparseMarker* marker, std::ostream *outStream);

#endif // SPARSEMARKER_H
