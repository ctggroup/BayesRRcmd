#ifndef LOGWRITER_H
#define LOGWRITER_H

#include <Eigen/Eigen>
#include <fstream>
#include <string>
#include "writer.h"

class LogWriter : public Writer 
{
public:
 
   void open();
 
};
#endif //LOGWRITER_H
