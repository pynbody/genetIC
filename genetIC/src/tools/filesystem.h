#ifndef IC_CHANGECWDWHILEINSCOPE_H
#define IC_CHANGECWDWHILEINSCOPE_H

#include <string>

namespace tools {
  class ChangeCwdWhileInScope {
  protected:
    std::string old;

  public:
    ChangeCwdWhileInScope(std::string newFolder);

    ~ChangeCwdWhileInScope();
  };

  std::string getDirectoryName(std::string full);
}

#endif
