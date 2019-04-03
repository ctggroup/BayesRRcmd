#include "densemarker.h"

#include "compression.h"

#include <fstream>

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


template<>
CompressedMarker compress(const DenseMarker *marker)
{
    CompressedMarker compressed;
    const auto maxCompressedOutputSize = maxCompressedDataSize<double>(marker->numInds);
    compressed.buffer.reset(new unsigned char[maxCompressedOutputSize]);
    compressed.size = compressData(*marker->Cx,
                                   compressed.buffer.get(),
                                   maxCompressedOutputSize);
    return compressed;
}

template<>
void write(const DenseMarker *marker, std::ostream *outStream)
{
    outStream->write(reinterpret_cast<char *>(marker->buffer.get()),
                     marker->Cx->size() * sizeof(double));
}


template<>
void MarkerBuilder<DenseMarker>::initialise(const unsigned int snp,
                                            const unsigned int numInds)
{
    reset(numInds);

    m_marker.reset(new DenseMarker);
    m_marker->i = snp;
    m_marker->numInds = numInds;

    const unsigned int colSize = numInds * sizeof (double);
    m_marker->buffer.reset(new unsigned char[colSize]);
    m_marker->Cx = std::make_shared<Map<VectorXd>>(reinterpret_cast<double *>(m_marker->buffer.get()),
                                                   numInds);
}

template<>
void MarkerBuilder<DenseMarker>::processAllele(unsigned int individual,
                                               unsigned int allele1,
                                               unsigned int allele2)
{
    if (allele1 == 0 && allele2 == 1) {  // missing genotype
        // Don't store a marker value like this as it requires floating point comparisons later
        // which are not done properly. Instead, store the index of the individual in a vector and simply
        // iterate over the collected indices. Also means iterating over far fewer elements which may
        // make a noticeable difference as this scales up.
        m_missingIndices.push_back(individual);
    } else {
        const auto value = allele1 + allele2;
        (*m_marker->Cx)[individual] = value;
        m_sum += value;
    }
}

template<>
void MarkerBuilder<DenseMarker>::endColumn()
{
    auto& snpData = *m_marker->Cx;

    const double mean = m_sum / (m_numInds - static_cast<double>(m_missingIndices.size()));
    std::for_each(m_missingIndices.cbegin(), m_missingIndices.cend(), [&](const unsigned int index) {
        snpData[index] = mean;
    });

    // Standardize genotypes
    snpData.array() -= snpData.mean();
    const auto sqn = snpData.squaredNorm();
    const auto sigma = 1.0 / (sqrt(sqn / (m_numInds - 1.0)));
    snpData.array() *= sigma;
}
