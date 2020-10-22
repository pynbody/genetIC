#ifndef _VECTORMATH_HPP_INCLUDED
#define _VECTORMATH_HPP_INCLUDED

#include <cassert>
#include <vector>
#include <complex>
#include <cmath>

namespace tools {
  namespace numerics {
    //! Multiplies vector a by constant b
    template<typename T, typename S>
    void operator*=(std::vector<T> &a, S b) {
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        a[i] *= b;
      }
    }

    //! Divides vector a by constant b
    template<typename T, typename S>
    void operator/=(std::vector<T> &a, S b) {
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        a[i] /= b;
      }
    }

    //! Multiplies each element of vector a by the corresponding element of vector b
    template<typename T>
    void operator*=(std::vector<T> &a, const std::vector<T> &b) {
      assert(a.size() == b.size());
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        a[i] *= b[i];
      }
    }

    //! Adds b to a, element-wise
    template<typename T>
    void operator+=(std::vector<T> &a, const std::vector<T> &b) {
      assert(a.size() == b.size());
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        a[i] += b[i];
      }
    }

    //! Divides each element of vector a by the corresponding element of vector b
    template<typename T>
    void operator/=(std::vector<T> &a, const std::vector<T> &b) {
      assert(a.size() == b.size());
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        a[i] /= b[i];
      }
    }

    /*! \brief Outputs the content of a vector to a stream
        \param output - stream to output vector to.
        \param a - vector to output.
    */
    template<typename T>
    std::ostream &operator<<(std::ostream &output, const std::vector<T> &a) {
      output << "[";
      for (auto a_element : a)
        output << a_element << ", ";
      output << "]";
      return output;
    }

    //! Returns a vector whose elements are the real parts of the elements of input.
    template<typename T>
    std::vector<T> real(const std::vector<std::complex<T>> &input) {
      std::vector<T> output;
      output.resize(input.size());
#pragma omp parallel for
      for (size_t i = 0; i < input.size(); ++i) {
        output[i] = std::real(input[i]);
      }
      return output;
    }

    //! For real inputs, returns a vector with the same elements as input
    template<typename T>
    std::vector<T> real(const std::vector<T> &input) {
      std::vector<T> output;
      output.resize(input.size());
#pragma omp parallel for
      for (size_t i = 0; i < input.size(); ++i) {
        output[i] = input[i];
      }
      return output;
    }

    //! Returns a vector whose elements are the absolute values of the elements of input
    template<typename T>
    std::vector<decltype(std::abs(T{}))> abs(const std::vector<T> &input) {
      std::vector<T> output;
      output.resize(input.size());
#pragma omp parallel for
      for (size_t i = 0; i < input.size(); ++i) {
        output[i] = std::abs(input[i]);
      }
      return output;
    }

    //! Returns a vector whose elements are the imaginary parts of the complex vector input.
    template<typename T>
    std::vector<T> imag(const std::vector<std::complex<T>> &input) {
      std::vector<T> output;
      output.resize(input.size());
#pragma omp parallel for
      for (size_t i = 0; i < input.size(); ++i) {
        output[i] = std::imag(input[i]);
      }
      return output;
    }

    //! Returns a vector whose elements are the products of the vectors a and b (which are not modified)
    template<typename T>
    std::vector<T> operator*(const std::vector<T> &a, const std::vector<T> &b) {
      std::vector<T> output;
      output.resize(a.size());
      assert(a.size() == b.size());
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        output[i] = a[i] * b[i];
      }
      return output;
    }

    //! Returns a vector whose elements are the ratio a_i/b_i for elements a_i and b_i of a and b respectively
    template<typename T>
    std::vector<T> operator/(const std::vector<T> &a, const std::vector<T> &b) {
      std::vector<T> output;
      output.resize(a.size());
      assert(a.size() == b.size());
#pragma omp parallel for
      for (size_t i = 0; i < a.size(); ++i) {
        output[i] = a[i] / b[i];
      }
      return output;
    }

    template<typename T>
    std::vector<T> log(const std::vector<T> &x) {
      std::vector<T> output;
      for(auto i = x.begin(); i!=x.end(); i++) {
        output.push_back(std::log(*i));
        logging::entry() << output.back() << " " << std::endl;
      }
      return output;
    }
  }
}
#endif
