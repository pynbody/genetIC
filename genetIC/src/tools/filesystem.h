#ifndef IC_CHANGECWDWHILEINSCOPE_H
#define IC_CHANGECWDWHILEINSCOPE_H

#include <string>

namespace tools {
  /*! \class ChangeCwdWhileInScope
      \brief Changes the current working directly. Usually used to run a different parameter file as an input mapper
  */
  class ChangeCwdWhileInScope {
  protected:
    std::string old; //!< String to store old directory, so that we can return to it.

  public:
    //! Stores the current directory, and then changes to the new one specified by newFolder
    ChangeCwdWhileInScope(std::string newFolder);

    //! Destructor - returns to the original working directory
    ~ChangeCwdWhileInScope();
  };

  //! Returns the directory of the given file
  std::string getDirectoryName(std::string full);
}

#endif
