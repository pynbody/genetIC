#ifndef _IO_INCLUDED
#define _IO_INCLUDED

//#include <cassert>
//#include <vector>
//#include <limits>

#include "src/simulation/particles/mapper/mapper.hpp"


namespace cosmology {
  template<typename FloatType>
  struct CosmologicalParameters;
}

#include "io/input.hpp"
#include "io/gadget.hpp"
#include "io/tipsy.hpp"
#include "io/grafic.hpp"

#include <iostream>
#include <string>
#include <stdexcept>

/*!
    \namespace io
    \brief Input output.

    Manage the different types of data outputs produced by the code. Possibilities are tipsy, grafic, gadget and so on
 */
namespace io {

  enum class OutputFormat {
    unknown = 1, gadget2 = 2, gadget3 = 3, tipsy = 4, grafic = 5
  };

  std::ostream &operator<<(std::ostream &outputStream, const OutputFormat &format) {
    switch (format) {
      case OutputFormat::unknown:
        outputStream << "unknown";
        break;
      case OutputFormat::gadget2:
        outputStream << "gadget2";
        break;
      case OutputFormat::gadget3:
        outputStream << "gadget3";
        break;
      case OutputFormat::tipsy:
        outputStream << "tipsy";
        break;
      case OutputFormat::grafic:
        outputStream << "grafic";
    }
    return outputStream;
  }

  std::istream &operator>>(std::istream &inputStream, OutputFormat &format) {
    std::string s;
    inputStream >> s;
    try {
      int i = std::stoi(s);
      format = static_cast<OutputFormat>(i);
    } catch (std::invalid_argument e) {
      if (s == "gadget2") {
        format = OutputFormat::gadget2;
      } else if (s == "gadget3") {
        format = OutputFormat::gadget3;
      } else if (s == "tipsy") {
        format = OutputFormat::tipsy;
      } else if (s == "grafic") {
        format = OutputFormat::grafic;
      } else {
        inputStream.setstate(std::ios::failbit);
      }
    }
    return inputStream;
  }




  /*! \namespace io::numpy
      \brief Deals with output to numpy formats.
  */

  /*! \namespace io::numpy::detail
      \brief Detailed implementation of numpy output/input.
  */

}


#endif
