//
// Created by Andrew Pontzen on 22/10/2020.
//

#include "logging.hpp"
#include <iostream>

namespace logging {
  namespace {
    unsigned numIndent = 0;
  }

  std::ostream & entry() {
    std::cerr << "";
    for(unsigned i=0; i<numIndent; i++)
      std::cerr << "| ";
    return std::cerr;
  }

  IndentWhileInScope::IndentWhileInScope() {
    numIndent+=1;
  }

  IndentWhileInScope::~IndentWhileInScope() {
    numIndent-=1;
  }


}