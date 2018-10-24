#ifndef IC_INPUT_HPP
#define IC_INPUT_HPP

#include "src/io.hpp"
#include <sstream>

namespace io {
  template<typename T>
  void getBuffer(std::vector<T> &store, std::string filename) {
    std::ifstream f(filename);
    if (!f.is_open())
      throw std::runtime_error("File " + filename + " not found");
    std::string line;
    while (!f.eof()) {
      T temp;
      if (f >> temp)
        store.push_back(temp);
      if (f.fail() && !f.eof())
        throw std::runtime_error("Error reading file " + filename);
    }
  }

  //! Returns the number of columns in the supplied data file.
  size_t getNumberOfColumns(std::string filename)
  {
    std::ifstream f(filename);
    std::stringstream ss;
    if(!f.is_open())
    {
        throw std::runtime_error("File " + filename + " not found");
    }
    //Extract the first line, and count the number of entries:
    std::string line;
    std::getline(f,line);
    ss << line;
    size_t nCols = 0;
    while(!ss.eof())
    {
        double temp;
        if(ss >> temp)
        {
            nCols++;
        }
        else
        {
            throw std::runtime_error("Error reading file " + filename + ". Could not determine number of columns.");
        }
    }
    return nCols;
  }

  template<typename T>
  void dumpBuffer(const std::vector<T> &store, std::string filename) {
    std::ofstream f(filename);
    if (!f.is_open())
      throw std::runtime_error("Can't open file " + filename);
    for (const T &item : store) {
      f << item << std::endl;
    }
  }
}

#endif //IC_INPUT_HPP
