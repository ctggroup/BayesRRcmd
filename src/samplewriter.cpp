#include "samplewriter.h"
#include "data.hpp"

#include <iostream>

using namespace Eigen;

void SampleWriter::open()
{
    std::cout << "Opening results file " << m_fileName << std::endl;
    m_outFile.open(m_fileName);

    m_outFile << "iteration," << "mu,";
    for (unsigned int i = 0; i < m_markerCount; ++i) {
        m_outFile << "beta[" << (i+1) << "],";
    }

    m_outFile << "sigmaE," << "sigmaG,";
    for (unsigned int i = 0; i < m_markerCount; ++i) {
        m_outFile << "comp[" << (i+1) << "],";
    }

    unsigned int i;
    for (i = 0; i < (m_individualCount - 1); ++i) {
        m_outFile << "epsilon[" << (i + 1) << "],";
    }

    m_outFile << "epsilon[" << (i + 1) << "]";
    m_outFile << std::endl;
    m_outFile.flush();
}

void SampleWriter::openGroups(int numberGroups)
{
    std::cout << "Opening results file " << m_fileName << std::endl;
    m_outFile.open(m_fileName);

    m_outFile << "iteration," << "mu,";
    for (unsigned int i = 0; i < m_markerCount; ++i) {
        m_outFile << "beta[" << (i+1) << "],";
    }

    m_outFile << "sigmaE,";

    for(unsigned int i = 0; i < numberGroups; ++i){
    	    m_outFile << "sigmaG[" << (i+1) << "],";
    	  }

    for (unsigned int i = 0; i < m_markerCount; ++i) {
        m_outFile << "comp[" << (i+1) << "],";
    }

    unsigned int i;
    for (i = 0; i < (m_individualCount - 1); ++i) {
        m_outFile << "epsilon[" << (i + 1) << "],";
    }

    m_outFile << "epsilon[" << (i + 1) << "]";
    m_outFile << std::endl;
    m_outFile.flush();
}
