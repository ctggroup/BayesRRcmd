#include "writer.h"

#include <iostream>
#include <Eigen/Eigen>

using namespace Eigen;


Writer::Writer() : m_commaInitFormat(StreamPrecision,DontAlignCols,", ", ", ", "", "", "", "")
{
}

Writer::~Writer()
{
  close();
}

void Writer::write(const Eigen::VectorXd &sample)
{
  m_outFile << sample.transpose().format(m_commaInitFormat) << std::endl;
  m_outFile.flush();
}

void Writer::close()
{
  m_outFile.close();
}
