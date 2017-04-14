#ifndef _UTILS_HPP_INCLUDED
#define _UTILS_HPP_INCLUDED

#include <vector>
/*!
    \namespace tools
    \brief Defines useful functions and tools used throughout the code

 */
namespace tools {
  // Argsort function from http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
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

  template<typename T>
  size_t getRatioAndAssertPositiveInteger(T p, T q, T tolerance = 1e-8) {
    assert(p > 0);
    assert(q > 0);
    T ratio = p / q;
    size_t rounded_ratio = size_t(round(ratio));
    assert(abs(T(rounded_ratio) - ratio) < tolerance);
    return rounded_ratio;
  }

  template<typename T>
  int getRatioAndAssertInteger(T p, T q, T tolerance = 1e-8) {
    T ratio = p / q;
    int rounded_ratio = int(round(ratio));
    assert(abs(T(rounded_ratio) - ratio) < tolerance);
    return rounded_ratio;
  }

  template<typename T>
  void set_zero(T &item) {
    item = 0;
  }

  template<typename T>
  void set_zero(std::complex<T> &item) {
    item.real(0);
    item.imag(0);
  }

  // helper function for dumps in mapper
  void indent(std::ostream &s, int level = 0) {
    for (int i = 0; i < level; i++) s << "| ";
    s << "+ ";

  }
}

#endif
