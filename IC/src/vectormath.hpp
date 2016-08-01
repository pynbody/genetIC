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


#endif