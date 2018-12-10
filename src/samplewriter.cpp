#include "samplewriter.h"

#include <iostream>
#include "options.hpp"

using namespace Eigen;

SampleWriter::SampleWriter()
    : m_commaInitFormat(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "")
{
}

SampleWriter::~SampleWriter()
{
    close();
}

void SampleWriter::open(std::string &header)
{
    std::cout << "Opening results file " << m_fileName << std::endl;
    m_outFile.open(m_fileName);
    m_outFile << header;
    m_outFile.flush();
}

void SampleWriter::write(const Eigen::VectorXd &sample)
{
    //std::cout << "Writing sample" << std::endl;
    m_outFile << sample.transpose().format(m_commaInitFormat) << std::endl;
    m_outFile.flush();
}

void SampleWriter::close()
{
    m_outFile.close();
}
