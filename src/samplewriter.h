#ifndef SAMPLEWRITER_H
#define SAMPLEWRITER_H

#include <Eigen/Eigen>
#include "writer.h"
#include <fstream>
#include <string>

class SampleWriter : public Writer
{
public:
    void setMarkerCount(unsigned int markerCount) { m_markerCount = markerCount; }
    unsigned int markerCount() const { return m_markerCount; }

    void setIndividualCount(unsigned int individualCount) { m_individualCount = individualCount; }
    unsigned int individualCount() const { return m_individualCount; }
    
    void open() override;
    void openGroups(int numberGroups);

    void setFixedCount(unsigned int fixedCount) { m_fixedCount = fixedCount; }
    unsigned int fixedCount() const { return m_fixedCount; }

    void open_bayesW();
    void open_bayesW_fixed();

private:
    unsigned int m_markerCount;
    unsigned int m_individualCount;
    unsigned int m_fixedCount;
};

#endif // SAMPLEWRITER_H
