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


/*!
    \namespace io
    \brief Input output.

    Manage the different types of data outputs produced by the code. Possibilities are tipsy, grafic, gadget and so on
 */
namespace io {

  enum class OutputFormat {
    gadget2 = 2, gadget3, tipsy, grafic
  };

  /*! \namespace io::numpy
      \brief Deals with output to numpy formats.
  */

  /*! \namespace io::numpy::detail
      \brief Detailed implementation of numpy output/input.
  */

}


#endif
