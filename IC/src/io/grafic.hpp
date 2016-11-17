//
// Created by Andrew Pontzen on 17/11/2016.
//

#include <fstream>
#include "../multilevelcontext.hpp"

template<typename T>
class GraficOutput {
protected:
  std::ofstream outputFile;
  const MultiLevelContextInformation<T> & context;
public:
  GraficOutput(std std::string & fname,
               const MultiLevelContextInformation<T> & context ) :
    outputFile(fname, std::ios::binary),
    context(context) {

  }

  void write() {
    context.forEachLevel([&]() {
      
    });
  }

};