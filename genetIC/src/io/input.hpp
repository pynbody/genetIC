//
// Created by Andrew Pontzen on 18/11/2016.
//

#ifndef IC_INPUT_HPP
#define IC_INPUT_HPP

#include "../io.hpp"

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
