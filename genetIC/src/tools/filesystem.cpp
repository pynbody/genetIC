#include "src/tools/filesystem.h"
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>
#include <iostream>
#include <stdexcept>

namespace tools {
    std::string getDirectoryName(std::string full) {
      char *fullCopy = new char[full.size() + 1];
      std::copy(full.begin(), full.end(), fullCopy);
      std::string rVal = dirname(fullCopy);
      delete[] fullCopy;
      return rVal;
    }


    ChangeCwdWhileInScope::ChangeCwdWhileInScope(std::string newFolder) {
      char buffer[1000];
      if (getcwd(buffer, 1000) == nullptr)
        throw std::runtime_error("Unable to get current working directory");
      old = buffer;
      chdir(newFolder.c_str());
    }

    ChangeCwdWhileInScope::~ChangeCwdWhileInScope() {
      chdir(old.c_str());
    }
}