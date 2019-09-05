#include "analysis.h"

#include <numeric>

Analysis::Analysis(const Data *data, const Options *opt)
    : m_data(data)
    , m_opt(opt)
{

}

Analysis::~Analysis() = default;

IndexEntry Analysis::indexEntry(unsigned int i) const
{
    if (!m_data)
        return {};

    return m_data->ppbedIndex[i];
}

bool Analysis::compressed() const
{
    return m_opt->compress;
}

unsigned char* Analysis::compressedData() const
{
    if (!m_data)
        return nullptr;

    return reinterpret_cast<unsigned char*>(m_data->ppBedMap);
}

std::string Analysis::preprocessedFile() const
{
    return ppFileForType(m_opt->preprocessDataType, m_opt->dataFile);
}

int Analysis::runGibbs(AnalysisGraph *analysis)
{
    std::vector<unsigned int> markerI(m_data->numSnps);
    std::iota(markerI.begin(), markerI.end(), 0);

    return runGibbs(analysis, std::move(markerI));
}

int Analysis::runGibbs(AnalysisGraph *analysis, const std::vector<unsigned int> &markers)
{
    const auto copy(markers);
    return runGibbs(analysis, std::move(copy));
}
