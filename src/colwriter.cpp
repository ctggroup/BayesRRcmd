#include "colwriter.h"
#include <iostream>

using namespace Eigen;


void ColWriter::open()
{
  std::cout << "Opening column log file" << m_fileName << std::endl;
  m_outFile.open(m_fileName);
  m_outFile << "marker,";
  m_outFile << "num_upd_microS,";
  m_outFile << "beta_upd_microS,";
  m_outFile << "eps_upd_microS";
  m_outFile << std::endl;
}
