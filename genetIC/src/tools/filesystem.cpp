#include <stdexcept>
#include "src/tools/filesystem.h"
#include <unistd.h>
#include <libgen.h>
#include <iostream>

namespace tools {
  //! Returns the directory of the given file
  std::string getDirectoryName(std::string full) {
    char *fullCopy = new char[full.size() + 1];
    std::copy(full.begin(), full.end(), fullCopy);
    std::string rVal = dirname(fullCopy);
    delete[] fullCopy;
    return rVal;
  }


  //! Stores the current directory, and then changes to the new one specified by newFolder
  ChangeCwdWhileInScope::ChangeCwdWhileInScope(std::string newFolder) {
    char buffer[1000];
    if (getcwd(buffer, 1000) == nullptr)
      throw std::runtime_error("Unable to get current working directory");
    old = buffer;
    chdir(newFolder.c_str());
  }

  //! Destructor - returns to the original working directory
  ChangeCwdWhileInScope::~ChangeCwdWhileInScope() {
    chdir(old.c_str());
  }
}
