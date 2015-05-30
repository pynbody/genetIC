#ifndef _ARGSORT_HPP_INCLUDED
#define _ARGSORT_HPP_INCLUDED

// Argsort function from http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

template <typename T>
vector<size_t> argsort(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
#endif
