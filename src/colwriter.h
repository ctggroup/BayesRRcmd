#ifndef COLWRITER_H
#define COLWRITER_H

#include "logwriter.h"
#include <iostream>
#include "writer.h"

class ColWriter : public Writer
{
 public:
  void open();
};

#endif //COLWRITER_H
