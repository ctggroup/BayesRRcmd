#ifndef SPARSEMARKER_H
#define SPARSEMARKER_H

#include "marker.h"

struct SparseMarker : public Marker
{
    double mean = 0;
    double sd = 0;
    double squrdZ= 0;
    double Zsum = 0;

    double numInds = 0;
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


struct EigenSparseMarker : public SparseMarker
{
    using UnitDataType = double;
    SparseVector<UnitDataType> Zg;

    const VectorXd *ones = nullptr;

    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;

protected:
    double dot(const VectorXd &epsilon) const override;
};


struct RaggedSparseMarker : public SparseMarker
{
    using IndexVector = std::vector<int>;
    // the indexes of elements of the bed matrix which are one for this column
    IndexVector Zones;
    // the indexes of elements of the bed matrix which are two for this column
    IndexVector Ztwos;

    // the indexes of elements of the bed matrix which are missing for this column
    IndexVector Zmissing;

    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;

protected:
    double dot(const VectorXd &epsilon) const override;
};

#endif // SPARSEMARKER_H
