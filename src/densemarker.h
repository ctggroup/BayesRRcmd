#ifndef DENSEMARKER_H
#define DENSEMARKER_H

#include "marker.h"
#include "markerbuilder.h"

struct DenseMarker : public Marker
{
    double component = 0;

    std::shared_ptr<unsigned char[]> buffer = nullptr;
    std::shared_ptr<Map<VectorXd>> Cx = nullptr;

    double computeNum(VectorXd &epsilon, const double beta_old) override;
    void updateEpsilon(VectorXd &epsilon,
                       const double beta_old,
                       const double beta) override;
};


template<>
CompressedMarker compress(const DenseMarker* marker);

template<>
void write(const DenseMarker* marker, std::ostream *outStream);

template<>
void MarkerBuilder<DenseMarker>::initialise(const unsigned int snp,
                                            const unsigned int numInds);

template<>
void MarkerBuilder<DenseMarker>::processAllele(unsigned int individual,
                                               unsigned int allele1,
                                               unsigned int allele2);
template<>
void MarkerBuilder<DenseMarker>::endColumn();

#endif // DENSEMARKER_H
