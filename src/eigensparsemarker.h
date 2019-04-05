#ifndef EIGENSPARSEMARKER_H
#define EIGENSPARSEMARKER_H

#include "sparsemarker.h"

struct EigenSparseMarker : public SparseMarker
{
    using UnitDataType = double;
    SparseVector<UnitDataType> Zg;

    const VectorXd *ones = nullptr;

    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;

    std::streamsize size() const override;
    void read(std::istream *inStream) override;
    void write(std::ostream *outStream) const override;

    bool isValid() const override;

protected:
    double dot(const VectorXd &epsilon) const override;
};

#endif // EIGENSPARSEMARKER_H
