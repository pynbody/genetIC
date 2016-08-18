//
// Created by Andrew Pontzen on 26/07/15.
//

#ifndef IC_FILTER_HPP
#define IC_FILTER_HPP

#include <math.h>
#include <stdexcept>



template<typename T>
class DivisionFilter;

template<typename T>
class Filter {
public:
  virtual ~Filter() { }

  virtual T operator()(T x) const {
    return 1.0;
  }

  template<typename S>
  auto compositeFunction(const S & other) const {
    using ret_type = decltype(other(T{}));
    auto func = [this, &other](T x) {
      return (*this)(x)*other(x);
    };
    return func;
  }

  virtual Filter<T> *clone() const {
    return new Filter<T>();
  }

  DivisionFilter<T> operator/(const Filter<T> &other) const {
    return DivisionFilter<T>(this, &other);
  }

  virtual void debugInfo(std::ostream &s) const {
    s << "Filter";
  }
};


template<typename T>
class DivisionFilter : public Filter <T> {
private:
  Filter<T> *pFirst, *pSecond;
public:
  DivisionFilter(const Filter<T> * first, const Filter<T> *second) : pFirst(first->clone()), pSecond(second->clone())
  {

  }

  virtual ~DivisionFilter() {
    delete pFirst;
    delete pSecond;
  }

  virtual T operator()(T x) const {
    return (*pFirst)(x)/(*pSecond)(x);
  }

  virtual void debugInfo(std::ostream &s) const override {
    s << (*pFirst) << "/" << (*pSecond);
  }
};



template<typename T>
class LowPassFermiFilter : public Filter<T> {
private:
  T kcut;
  T temperature;

public:
  LowPassFermiFilter(T kcut_in) : kcut(kcut_in) {
    // The following arbitrary choice for the temperature (i.e.
    // sharpness of the cut-off) has been made by numerical experiments
    // trading off the various pros and cons of a sharp cut discussed
    // in the notes.
    temperature = kcut / 10;
  };

  virtual Filter<T> *clone() const {
    return new LowPassFermiFilter<T>(*this);
  }

  T operator()(T k) const override {
    return 1. / (1. + exp((k - kcut) / temperature));
  }

  virtual void debugInfo(std::ostream &s) const override {
    s << "LowPassFermiFilter(" << kcut << ")";
  }
};

template<typename T>
class ComplementaryCovarianceFilterAdaptor : public Filter<T> {
private:
  Filter<T> *pUnderlying;
public:
  ComplementaryCovarianceFilterAdaptor(const Filter<T> &original) : pUnderlying(original.clone()) { };

  virtual ~ComplementaryCovarianceFilterAdaptor() {
    delete pUnderlying;
  }

  T operator()(T k) const override {
    T denFilter = (*pUnderlying)(k);
    return sqrt(1.-denFilter*denFilter);
  }

  virtual Filter<T> *clone() const {
    return new ComplementaryCovarianceFilterAdaptor<T>(*this);
  }

  virtual void debugInfo(std::ostream &s) const override {
    s << "ComplementaryCovarianceFilterAdaptor(" << (*pUnderlying) << ")";
  }
};

template<typename T>
class ComplementaryFilterAdaptor : public Filter<T> {
private:
  Filter<T> *pUnderlying;
public:
  ComplementaryFilterAdaptor(const Filter<T> &original) : pUnderlying(original.clone()) { };

  virtual ~ComplementaryFilterAdaptor() {
    delete pUnderlying;
  }

  T operator()(T k) const override {
    return 1. - (*pUnderlying)(k);
  }

  virtual Filter<T> *clone() const {
    return new ComplementaryFilterAdaptor<T>(*this);
  }

  virtual void debugInfo(std::ostream &s) const override {
    s << "ComplementaryFilterAdaptor(" << (*pUnderlying) << ")";
  }
};


template<typename T>
class FilterFamily {
protected:
  std::vector<Filter<T> *> filters;

  FilterFamily() {

  }

public:

  FilterFamily(size_t maxLevels) {
    for(size_t i=0; i<maxLevels; ++i)
      filters.push_back(new Filter<T>());
  }

  virtual ~FilterFamily() {
    /*
    for(auto f: filters)
      delete f;*/
  };

  virtual const Filter<T> &getFilterOnLevel(int level) const {
    return *(filters[level]);
  }

  virtual size_t getMaxLevel() const {
    return filters.size();
  }

  virtual void debugInfo(std::ostream &s) const  {
    s << "FilterFamily(";
    for(size_t i=0; i<filters.size(); ++i) {
      s << *filters[i];
      if(i+1<filters.size()) s << ", ";
    }
    s << ")";

  }
};

template<typename T>
class TwoLevelFilterFamily : public FilterFamily<T> {
protected:

public:

  TwoLevelFilterFamily(T k_cut)
  {
    this->filters.push_back(new LowPassFermiFilter<T>(k_cut));
    this->filters.push_back(new ComplementaryCovarianceFilterAdaptor<T>(LowPassFermiFilter<T>(k_cut)));
  }

};


template <typename T>
class ResidualFilterFamily : public FilterFamily<T> {

public:
  ResidualFilterFamily(const FilterFamily<T> &underlying)
  {
    for(size_t i=0; i<underlying.getMaxLevel(); i++) {
      this->filters.push_back(new ComplementaryFilterAdaptor<T>(underlying.getFilterOnLevel(i)));
    }
  }

};



template<typename T>
class TwoLevelDependentFilterFamily : public FilterFamily<T> {
public:

  TwoLevelDependentFilterFamily(T k_cut)
  {
    this->filters.push_back(new LowPassFermiFilter<T>(k_cut));
    this->filters.push_back(new ComplementaryFilterAdaptor<T>(LowPassFermiFilter<T>(k_cut)));
  }

};

template<typename T>
class TwoLevelRecombinedFilterFamily : public FilterFamily<T> {
public:

  TwoLevelRecombinedFilterFamily(T k_cut)
  {
    this->filters.push_back(new ComplementaryFilterAdaptor<T>(Filter<T>()));
    this->filters.push_back(new Filter<T>());
  }

};


template<typename T>
std::ostream & operator<<(std::ostream &s, const Filter<T> &f) {
  f.debugInfo(s);
  return s;
}

template<typename T>
std::ostream & operator<<(std::ostream &s, const FilterFamily<T> &f) {
  f.debugInfo(s);
  return s;
}

template<typename T>
void tabulateFilter(Filter<T> *pF) {
  std::cerr << "tabulateFilter:" << std::endl;
  for (T k = 0.1; k < 3.0; k += 0.1) {
    std::cerr << " " << k << " " << (*pF)(k) << std::endl;
  }
}


#endif //IC_FILTER_HPP
