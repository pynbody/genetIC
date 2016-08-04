#ifndef _SPARSE_HPP_INCLUDED
#define _SPARSE_HPP_INCLUDED

#include <cassert>
#include <vector>

template<typename T, typename S>
void operator*=(std::vector<T> &a, S b) {
#pragma omp parallel for
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] *= b;
  }
}

template<typename T, typename S>
void operator/=(std::vector<T> &a, S b) {
#pragma omp parallel for
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] /= b;
  }
}

template<typename T>
void operator*=(std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size()==b.size());
#pragma omp parallel for
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] *= b[i];
  }
}

template<typename T>
void operator/=(std::vector<T> &a, const std::vector<T> &b) {
  assert(a.size()==b.size());
#pragma omp parallel for
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] /= b[i];
  }
}

template<typename T>
std::ostream & operator<<(std::ostream & output, const std::vector<T> & a) {
  for(auto a_element : a)
    output << a_element << " ";
  return output;
}

template<typename T>
std::vector<T> real(std::vector<std::complex<T>> input) {
  std::vector<T> output;
  output.resize(input.size());
#pragma omp parallel for
  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = std::real(input[i]);
  }
  return output;
}

template<typename T>
std::vector<T> imag(std::vector<std::complex<T>> input) {
  std::vector<T> output;
  output.resize(input.size());
#pragma omp parallel for
  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = std::imag(input[i]);
  }
  return output;
}

#endif