#include "logwriter.h"
#include <iostream>

using namespace Eigen;

void LogWriter::open()
{
    std::cout << "Opening iteration log file " << m_fileName << std::endl;
    m_outFile.open(m_fileName);
    m_outFile << "iter,";
    m_outFile << "m_0,";
    m_outFile << "sigmaG,";
    m_outFile << "sigmaE,";
    m_outFile << "mu,";
    m_outFile << "mu_upd_microS,";
    m_outFile << "sigmaG_upd_microS,";
    m_outFile << "sigmaE_upd_microS,";
    m_outFile << "flowgraph_milliS,";
    m_outFile << "iter_d_milliS";
    
    m_outFile << std::endl;
    m_outFile.flush();

}
