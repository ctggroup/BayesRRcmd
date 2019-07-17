#ifndef WRITER_H
#define WRITER_H

#include <Eigen/Eigen>
#include <fstream>
#include <string>

class Writer
{
 public:
  Writer();
  ~Writer();

  void setFileName(const std::string &fileName) {m_fileName = fileName;}
  std::string fileName() const { return m_fileName;}

  virtual void open()=0;
  void write(const Eigen::VectorXd &message);
  void close();

 protected:
  std::string m_fileName;
  std::ofstream m_outFile;
  Eigen::IOFormat m_commaInitFormat;
 
};

#endif //WRITER_H
