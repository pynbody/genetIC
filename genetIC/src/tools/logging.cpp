//
// Created by Andrew Pontzen on 22/10/2020.
//

#include "logging.hpp"
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <cstdio>

namespace logging {
  namespace {
    unsigned numIndent = 0;

    const char* RESET = "\x1B[0m";
    const char* RED = "\x1B[31m";
    const char* YELLOW = "\x1B[33m";
    const char* BLUE = "\x1B[34m";


  }

  std::ostream & entry(level lev) {
    const bool isTty = isatty(fileno(stderr));

    if(isTty)
      std::cerr << RESET;

    std::time_t t = std::time(nullptr);

    std::cerr << std::put_time(std::localtime(&t), " %F %T   ");

    for(unsigned i=0; i<numIndent; i++)
      std::cerr << "| ";

    if(isTty && lev == level::warning)
      std::cerr << YELLOW;

    return std::cerr;
  }

  IndentWhileInScope::IndentWhileInScope() {
    numIndent+=1;
  }

  IndentWhileInScope::~IndentWhileInScope() {
    numIndent-=1;
  }


}