//
// Created by Andrew Pontzen on 22/10/2020.
//

#ifndef IC_LOGGING_HPP
#define IC_LOGGING_HPP

#include <ostream>

namespace logging {
  std::ostream & entry();

  class IndentWhileInScope {
  public:
    IndentWhileInScope();
    ~IndentWhileInScope();
  };

}
#endif //IC_LOGGING_HPP
