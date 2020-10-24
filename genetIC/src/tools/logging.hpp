//
// Created by Andrew Pontzen on 22/10/2020.
//

#ifndef IC_LOGGING_HPP
#define IC_LOGGING_HPP

#include <ostream>

namespace logging {
  enum level {
    warning = 0, info = 1, debug = 2
  };

  std::ostream & entry(level lev = info);

  class IndentWhileInScope {
  public:
    IndentWhileInScope();
    ~IndentWhileInScope();
  };

}
#endif //IC_LOGGING_HPP
