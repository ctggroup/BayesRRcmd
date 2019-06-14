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

private:
    unsigned int m_markerCount;
    unsigned int m_individualCount;
};

#endif // SAMPLEWRITER_H
