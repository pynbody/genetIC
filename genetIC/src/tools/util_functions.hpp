#ifndef _UTILS_HPP_INCLUDED
#define _UTILS_HPP_INCLUDED

#include <vector>
#include <cmath>
#include <stdexcept>
/*!
    \namespace tools
    \brief Defines useful functions and tools used throughout the code

 */
namespace tools {
  //! Argsort function from http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
  template<typename T>
  std::vector<size_t> argsort(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
  }

  /*! \brief Returns the integer ratio of two variables, throwing a run-time error if the ratio is not sufficiently close to an integer
      \param p - numerator
      \param q - denominator
      \param tolerance - how close to an integer we are willing to tolerate.
  */
  template<typename T>
  int getRatioAndAssertInteger(T p, T q, T tolerance = 1e-8) {
    T ratio = p / q;
    int rounded_ratio = int(round(ratio));
    if (!(std::abs(T(rounded_ratio) - ratio) < tolerance)){
      throw std::runtime_error("The ratio is not an integer within tolerance");
    }
    return rounded_ratio;
  }

  //! As getRatioAndAssertInteger, but additionally throws an error if the arguments are not positive.
  template<typename T>
  size_t getRatioAndAssertPositiveInteger(T p, T q, T tolerance = 1e-8) {
    assert(p > 0);
    assert(q > 0);
    return (size_t) getRatioAndAssertInteger(p, q, tolerance);
  }

  //! Solves n^x = p for integers x
  int findPowerOf(size_t n, size_t p){
    return getRatioAndAssertInteger(log(p), log(n));
  }

  //! Sets a given variable to zero (for real numbers)
  template<typename T>
  void set_zero(T &item) {
    item = 0;
  }

  //! Sets a given complex number to zero
  template<typename T>
  void set_zero(std::complex<T> &item) {
    item.real(0);
    item.imag(0);
  }

  //! Helper function for dumps in mapper
  void indent(std::ostream &s, int level = 0) {
    for (int i = 0; i < level; i++) s << "| ";
    s << "+ ";

  }

  //! Creates a linearly spaced vector with n elments between start and end.
  template<typename T>
  std::vector<T> linspace(T start, T end, int n) {
    std::vector<T> array;
    if ((n == 0) || (n == 1) || (start == end))
      array.push_back(end);
    else if (n > 1) {
      double step = (end - start) / (n);
      int count = 0;
      while (count < n + 1) {
        T elem = start + count * step;
        array.push_back(elem);
        ++count;
      }
    }
    return array;
  }

  //! Sorts a vector and removes any duplicated entries (used for sorting lists of flagged cells)
  template<typename T>
  void sortAndEraseDuplicate( std::vector<T> & vector){
    std::sort(vector.begin(), vector.end());
    vector.erase(std::unique(vector.begin(), vector.end()), vector.end());
  }


}

#endif
