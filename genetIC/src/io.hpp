#ifndef _IO_INCLUDED
#define _IO_INCLUDED

#include <cassert>
#include <vector>
#include <limits>

#include "src/simulation/particles/mapper/mapper.hpp"


#ifdef HAVE_HDF5
hid_t hdf_float = H5Tcopy (H5T_NATIVE_DOUBLE);
hid_t hdf_double = H5Tcopy (H5T_NATIVE_DOUBLE);
#endif


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

  enum class OutputFormat { gadget2=2, gadget3, tipsy, grafic };

}





#endif // _IO_INCLUDED
