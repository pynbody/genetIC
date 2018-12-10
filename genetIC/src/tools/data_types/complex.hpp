#ifndef IC_COMPLEX_HPP_H
#define IC_COMPLEX_HPP_H
/*!
    \namespace tools::datatypes
    \brief Provides definition for complex quantities and float precision

 */
namespace tools {
  namespace datatypes {

    //! Class template designed to ensure strip complex will return T if T is real
    template<typename T>
    class strip_complex_s {
    public:
      using type = T;
    };

    //! Class template designed to strip complexity if template parameter is complex
    template<typename T>
    class strip_complex_s<std::complex<T>> {
    public:
      using type = T;
    };

    //! Similar to strip_complex_s, but won't work with non-complex types (no type declared)
    template<typename T>
    class underlying_if_complex_s {

    };

    //! Similar to strip_complex_s, but won't work with non-complex types
    template<typename T>
    class underlying_if_complex_s<std::complex<T>> {
    public:
      using type = T;
    };


    //! Typename to easily remove complexity from a class
    template<typename T>
    using strip_complex = typename strip_complex_s<T>::type;

    //! As strip_complex, but only works to return the underlying type if the argument is actually complex
    template<typename T>
    using underlying_if_complex = typename underlying_if_complex_s<T>::type;

    //! Removes complexity if it exists, and then re-applies it, to guarantee a complex result.
    template<typename T>
    using ensure_complex = std::complex<strip_complex<T>>;

    //! Returns the real part of T (this specialisation is only reached if T is non-complex, so just returns the input)
    template<typename T>
    T real_part_if_complex(const T &input) {
      return input; // no-op for real types
    }

    //! Returns the real part of T (this specialisation only works for complex inputs, so .real() is always understood
    template<typename T>
    T real_part_if_complex(const std::complex<T> &input) {
      return input.real();
    }
  }
}

#endif //IC_COMPLEX_HPP_H
