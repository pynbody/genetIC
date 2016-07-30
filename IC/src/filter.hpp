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

  virtual Filter<T> *clone() const {
    return new Filter<T>();
  }

  DivisionFilter<T> operator/(const Filter<T> &other) const {
    return DivisionFilter<T>(this, &other);
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
};


template<typename T>
void tabulateFilter(Filter<T> *pF) {
  std::cerr << "tabulateFilter:" << std::endl;
  for (T k = 0.1; k < 3.0; k += 0.1) {
    std::cerr << " " << k << " " << (*pF)(k) << std::endl;
  }
}

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
};

template<typename T>
class ComplementaryCovarianceFilterAdaptor : public Filter<T> {
private:
  Filter<T> *pUnderlying;
public:
  ComplementaryCovarianceFilterAdaptor(Filter<T> *pOriginal) : pUnderlying(pOriginal) { };

  T operator()(T k) const override {
    T denFilter = (*pUnderlying)(k);
    return sqrt(1.-denFilter*denFilter);
  }

  virtual Filter<T> *clone() const {
    return new ComplementaryCovarianceFilterAdaptor<T>(*this);
  }
};

template<typename T>
class ComplementaryFilterAdaptor : public Filter<T> {
private:
  Filter<T> *pUnderlying;
public:
  ComplementaryFilterAdaptor(Filter<T> *pOriginal) : pUnderlying(pOriginal) { };

  T operator()(T k) const override {
    return 1. - (*pUnderlying)(k);
  }

  virtual Filter<T> *clone() const {
    return new ComplementaryFilterAdaptor<T>(*this);
  }
};


template<typename T>
class FilterFamily {
  Filter<T> filter;
public:
  virtual ~FilterFamily() { };

  virtual Filter<T> &getFilterOnLevel(int level) {
    return filter;
  }

  virtual Filter<T> &getResidualFilterOnLevel(int level) {
    throw std::runtime_error("FilterFamily cannot produce filter for requested level");
  }
};

template<typename T>
class TwoLevelFilterFamily : public FilterFamily<T> {
protected:
  Filter<T> *pDenFilterLow;
  Filter<T> *pDenFilterHigh;
  Filter<T> *pCovFilterHighPartOfLowRes;
  Filter<T> *pCovFilterHigh;

  TwoLevelFilterFamily() {

  }

public:



  TwoLevelFilterFamily(T k_cut) {
    pDenFilterLow = new LowPassFermiFilter<T>(k_cut);
    pDenFilterHigh = new ComplementaryCovarianceFilterAdaptor<T>(pDenFilterLow);
    pCovFilterHighPartOfLowRes = new ComplementaryFilterAdaptor<T>(pDenFilterLow);
  }

  ~TwoLevelFilterFamily() {
    delete pDenFilterLow;
    delete pDenFilterHigh;
    delete pCovFilterHighPartOfLowRes;
  }

  Filter<T> &getFilterOnLevel(int level) override {
    switch (level) {
      case 0:
        return *pDenFilterLow;
      case 1:
        return *pDenFilterHigh;
      default:
        throw std::runtime_error("TwoLevelFilterFamily cannot produce filter for requested level");
    }
  }


  Filter<T> &getResidualFilterOnLevel(int level) override {
    switch (level) {
      case 0:
        return *pCovFilterHighPartOfLowRes;
      default:
        throw std::runtime_error("TwoLevelFilterFamily cannot produce filter for requested level");
    }
  }
};




template<typename T>
class TwoLevelDependentFilterFamily : public TwoLevelFilterFamily<T> {
public:

  TwoLevelDependentFilterFamily(T k_cut) {
    this->pDenFilterLow = new LowPassFermiFilter<T>(k_cut);
    this->pDenFilterHigh = new ComplementaryFilterAdaptor<T>(this->pDenFilterLow);
    this->pCovFilterHighPartOfLowRes = new ComplementaryFilterAdaptor<T>(this->pDenFilterLow);
  }

};

#endif //IC_FILTER_HPP
