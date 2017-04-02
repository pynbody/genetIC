#ifndef IC_COMPLEX_HPP_H
#define IC_COMPLEX_HPP_H
/*!
    \namespace tools::datatypes
    \brief Provides definition for complex quantities and float precision

 */
namespace  tools {
    namespace datatypes {
        template<typename T>
        class strip_complex_s {
        public:
            using type = T;
        };

        template<typename T>
        class strip_complex_s<std::complex<T>> {
        public:
            using type = T;
        };

        template<typename T>
        class underlying_if_complex_s {

        };

        template<typename T>
        class underlying_if_complex_s<std::complex<T>> {
        public:
            using type = T;
        };


        template<typename T>
        using strip_complex = typename strip_complex_s<T>::type;

        template<typename T>
        using underlying_if_complex = typename underlying_if_complex_s<T>::type;

        template<typename T>
        using ensure_complex = std::complex<strip_complex<T>>;

        template<typename T>
        T real_part_if_complex(const T &input) {
            return input; // no-op for real types
        }

        template<typename T>
        T real_part_if_complex(const std::complex<T> &input) {
            return input.real();
        }
    }
}

#endif //IC_COMPLEX_HPP_H
