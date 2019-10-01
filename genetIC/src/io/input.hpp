#ifndef IC_INPUT_HPP
#define IC_INPUT_HPP

#include "src/io.hpp"
#include <sstream>

namespace io {

/*! \brief Reads data from a specified stream into a vector, until the end of file is reached

    \param store - vector in which to store the data.
    \param f - stream from which to extract data.
    \param filename - name of file that data is being extracted from, used as information in case of error
*/
  template<typename T>
  void getFromExistingStream(std::vector<T> &store, std::istream &f, std::string filename) {
    while (!f.eof()) {
      T temp;
      if (f >> temp)
        store.push_back(temp);
      if (f.fail() && !f.eof())
        throw std::runtime_error("Error reading file " + filename);
    }
  }

//! \brief Read numerical data from a file into a vector.
/*!
    Data must be formatted in columns separated by spaces, and rows
    separated by new lines. Data must also not contain any column
    or row labels - only numerical data.
*/
  template<typename T>
  void getBuffer(std::vector<T> &store, std::string filename) {
    std::ifstream f(filename);
    if (!f.is_open())
      throw std::runtime_error("File " + filename + " not found");

    getFromExistingStream(store, f, filename);
  }

  /*! \brief Process the data file, removing any potential column headers.

    Some versions of CAMB output column headers indicating what data each column refers to.
    These cannot be interpreted by getBuffer, so we need to remove them before processing.

  */
  template<typename T>
  void getBufferIgnoringColumnHeaders(std::vector<T> &store, std::string filename) {
    std::ifstream f(filename);
    std::stringstream ss;
    if (!f.is_open())
      throw std::runtime_error("File " + filename + " not found");

    std::string line;
    // Check the first line, to see if it contains valid numbers:
    std::getline(f, line);
    ss << line;
    T temp;
    if (ss >> temp) {
      // Can interpret first line as numbers, so continue as normal.
      // Finish the rest of this line first, then do the rest:
      store.push_back(temp);
      getFromExistingStream(store, ss, filename + " on first line.");
    }
    // Either way, we carry on with the rest of the file. From now on, errors should
    // be interpreted as corrupt data:
    getFromExistingStream(store, f, filename);
  }

  //! \brief Returns the number of columns in the supplied data file.
  //! This is used to distinguish between old and new CAMB formats.
  size_t getNumberOfColumns(std::string filename) {
    std::ifstream f(filename);
    std::stringstream ss;
    if (!f.is_open()) {
      throw std::runtime_error("File " + filename + " not found");
    }
    // Extract the second line, and count the number of entries...
    // If there is no second line, then the CAMB file can't be used for interpolation anyway
    std::string line;
    std::getline(f, line);
    std::getline(f, line);

    ss << line;
    size_t nCols = 0;
    while (!ss.eof()) {
      double temp;
      if (ss >> temp) {
        nCols++;
      } else {
        throw std::runtime_error("Error reading file " + filename + ". Could not determine number of columns.");
      }
    }
    return nCols;
  }

  //! \brief Output the numerical contents of the vector store to filename
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
