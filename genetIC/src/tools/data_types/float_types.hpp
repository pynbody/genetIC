#ifndef __FLOAT_TYPES
#define __FLOAT_TYPES

namespace tools {
  namespace datatypes {

    /*! \struct floatinfo
        \brief Base template defining information about floating point types that can be used by the code.

        This template will only be reached by types that aren't double or single precision. They are assumed to
        be non-double.
    */
    template<typename FloatType>
    struct floatinfo {
      static constexpr char const *name = "unknown";
      static constexpr int doubleprecision = 0;
    };

    //! Template specialisation for double precision floats
    template<>
    struct floatinfo<double> {
      static constexpr char const *name = "doub";
      static constexpr int doubleprecision = 1;
    };

    //! Template specialisation for single precision floats
    template<>
    struct floatinfo<float> {
      static constexpr char const *name = "float";
      static constexpr int doubleprecision = 0;
    };
  }
}

#endif
