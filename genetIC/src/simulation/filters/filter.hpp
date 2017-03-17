//
// Created by Andrew Pontzen on 26/07/15.
//

#ifndef IC_FILTER_HPP
#define IC_FILTER_HPP

#include <math.h>
#include <stdexcept>
/*!
    \namespace filters
    \brief Define the filters used to separate low and high k modes in the box

 */
namespace filters {

  template<typename T>
  class DivisionFilter;

  template<typename T>
  class ProductFilter;

  template<typename T>
  class Filter : public std::enable_shared_from_this<Filter<T>> {
  public:
    typedef T ReturnType;

    virtual ~Filter() {}

    Filter() {}

    Filter(T k_cut) {
      // provided for compatibility with filters that do take a wavenumber cut
    }

    virtual T operator()(T x) const {
      return 1.0;
    }

    template<typename S>
    auto compositeFunction(const S &other) const {
      auto func = [this, &other](T x) {
        return (*this)(x) * other(x);
      };
      return func;
    }

    virtual std::shared_ptr<Filter<T>> clone() const {
      return std::make_shared<Filter<T>>();
    }

    DivisionFilter<T> operator/(const Filter<T> &other) const {
      return DivisionFilter<T>(*this, other);
    }

    ProductFilter<T> operator*(const Filter<T> &other) const {
      return ProductFilter<T>(*this, other);
    }

    virtual void debugInfo(std::ostream &s) const {
      s << "Filter";
    }
  };

  template<typename T>
  class NullFilter : public Filter<T> {
  public:
    NullFilter() {


    }

    NullFilter(T k_cut) {
      // provided for compatibility with filters that do take a wavenumber cut
    }

    T operator()(T x) const override {
      return 0.0;
    }

    std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<NullFilter<T>>();
    }

    virtual void debugInfo(std::ostream &s) const {
      s << "NullFilter";
    }

  };


  template<typename T>
  class DivisionFilter : public Filter<T> {
  protected:
    std::shared_ptr<Filter<T>> pFirst, pSecond;
  public:
    DivisionFilter(const Filter<T> &first, const Filter<T> &second) : pFirst(first.clone()), pSecond(second.clone()) {

    }

    virtual T operator()(T x) const {
      return (*pFirst)(x) / (*pSecond)(x);
    }

    virtual std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<DivisionFilter<T>>(*pFirst, *pSecond);
    }

    virtual void debugInfo(std::ostream &s) const override {
      s << (*pFirst) << "/" << (*pSecond);
    }
  };

  template<typename T>
  class ProductFilter : public DivisionFilter<T> {
  protected:
    using DivisionFilter<T>::pFirst;
    using DivisionFilter<T>::pSecond;

  public:

    ProductFilter(const Filter<T> &first, const Filter<T> &second) : DivisionFilter<T>(first, second) {

    }

    virtual std::shared_ptr<Filter<T>> clone() const override {
      // simplification/optimization:
      if (typeid(*pFirst) == typeid(Filter<T>))
        return pSecond->clone();
      else if (typeid(*pSecond) == typeid(Filter<T>))
        return pFirst->clone();
      else
        return std::make_shared<ProductFilter<T>>(*pFirst, *pSecond);
    }

    virtual T operator()(T x) const {
      return (*pFirst)(x) * (*pSecond)(x);
    }

    virtual void debugInfo(std::ostream &s) const override {
      s << (*pFirst) << "*" << (*pSecond);
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

    virtual std::shared_ptr<Filter<T>> clone() const {
      return std::make_shared<LowPassFermiFilter<T>>(*this);
    }

    T operator()(T k) const override {
      return 1. / (1. + exp((k - kcut) / temperature));
    }

    virtual void debugInfo(std::ostream &s) const override {
      s << "LowPassFermiFilter(" << kcut << ")";
    }
  };

  template<typename UnderlyingType, typename T=typename UnderlyingType::ReturnType>
  class ComplementaryCovarianceFilterAdaptor : public Filter<T> {
  private:

    std::shared_ptr<Filter<T>> pUnderlying;
  public:

    template<typename... Args>
    ComplementaryCovarianceFilterAdaptor(Args &&... args) {
      pUnderlying = std::make_shared<UnderlyingType>(std::forward<Args>(args)...);
    };


    T operator()(T k) const override {
      T denFilter = (*pUnderlying)(k);
      return sqrt(1. - denFilter * denFilter);
    }

    virtual std::shared_ptr<Filter<T>> clone() const {
      return std::make_shared<ComplementaryCovarianceFilterAdaptor<UnderlyingType>>(*this);
    }

    virtual void debugInfo(std::ostream &s) const override {
      s << "ComplementaryCovarianceFilterAdaptor(" << (*pUnderlying) << ")";
    }
  };

  template<typename UnderlyingType, typename T=typename UnderlyingType::ReturnType>
  class ComplementaryFilterAdaptor : public Filter<T> {
  private:
    std::shared_ptr<Filter<T>> pUnderlying;
  public:

    template<typename... Args>
    ComplementaryFilterAdaptor(Args &&... args) {
      pUnderlying = std::make_shared<UnderlyingType>(std::forward<Args>(args)...);
    };

    ComplementaryFilterAdaptor(const Filter<T> &pOriginal) {
      assert(
        typeid(pOriginal) != typeid(ComplementaryFilterAdaptor<UnderlyingType, T>)); // this is NOT a copy constructor!
      pUnderlying = pOriginal.clone();
    }


    T operator()(T k) const override {
      return 1. - (*pUnderlying)(k);
    }

    virtual std::shared_ptr<Filter<T>> clone() const {
      return std::make_shared<ComplementaryFilterAdaptor<UnderlyingType>>(*this);
    }

    virtual void debugInfo(std::ostream &s) const override {
      s << "ComplementaryFilterAdaptor(" << (*pUnderlying) << ")";
    }
  };


  template<typename T>
  class ResidualFilterFamily;

  template<typename T>
  class FilterFamily {
  protected:
    friend class ResidualFilterFamily<T>;

    std::vector<std::shared_ptr<Filter<T>>> filters;
    std::vector<std::shared_ptr<Filter<T>>> complementFilters;
    std::vector<std::shared_ptr<Filter<T>>> hpFilters;
    std::vector<std::shared_ptr<Filter<T>>> lpFilters;

    FilterFamily() {

    }

  public:

    FilterFamily(size_t maxLevels) {
      for (size_t i = 0; i < maxLevels; ++i)
        filters.push_back(std::make_shared<Filter<T>>());
    }


    virtual const Filter<T> &getFilterOnLevel(int level) const {
      return *(filters[level]);
    }

    virtual const Filter<T> &getHighPassFilterOnLevel(int level) const {
      return *(hpFilters[level]);
    }

    virtual const Filter<T> &getLowPassFilterOnLevel(int level) const {
      return *(lpFilters[level]);
    }

    virtual size_t getMaxLevel() const {
      return filters.size();
    }

    virtual void debugInfo(std::ostream &s) const {
      s << "FilterFamily(";
      for (size_t i = 0; i < filters.size(); ++i) {
        s << *filters[i];
        if (i + 1 < filters.size()) s << ", ";
      }
      s << ")";

    }

    virtual void addLevel(T k_cut) {
      throw std::runtime_error("Don't know how to add a level to this type of filter family");
    }
  };

  template<typename LowFilterType, typename HighFilterType>
  class GenericMultiLevelFilterFamily : public FilterFamily<typename LowFilterType::ReturnType> {
  protected:
    using T = typename LowFilterType::ReturnType;
    static_assert(std::is_same<typename LowFilterType::ReturnType, typename HighFilterType::ReturnType>::value,
                  "Filters must use same floating point type");
  public:

    GenericMultiLevelFilterFamily() {
      this->filters.push_back(std::make_shared<Filter<T>>());
      this->lpFilters.push_back(std::make_shared<Filter<T>>());
      this->hpFilters.push_back(std::make_shared<Filter<T>>());
    }

    void addLevel(T k_cut) override {
      std::shared_ptr<Filter<T>> filt = this->filters.back();
      this->filters.pop_back();
      this->lpFilters.pop_back();

      // the low-pass filtered version of the field on what used to be the finest level:
      this->filters.push_back(((*filt) * LowFilterType(k_cut)).clone());

      this->complementFilters.push_back(((*filt) * HighFilterType(k_cut)).clone());

      // the low-pass filter to be used when copying from what used to be the finest level to the new finest level:
      this->lpFilters.push_back(LowFilterType(k_cut).clone());

      // the new finest level filter:
      this->filters.push_back(((*filt) * HighFilterType(k_cut)).clone());

      // the high-pass filter for removing unwanted low-k information from the new fine level:
      this->hpFilters.push_back(HighFilterType(k_cut).clone());

      // if this ends up being the last level, there is no low-pass filter ever applied
      this->lpFilters.push_back(Filter<T>().clone());

    }

  };

  template<typename T>
  using MultiLevelFilterFamily = GenericMultiLevelFilterFamily<LowPassFermiFilter<T>,
    ComplementaryCovarianceFilterAdaptor<LowPassFermiFilter<T>>>;

  template<typename T>
  using MultiLevelDependentFilterFamily = GenericMultiLevelFilterFamily<LowPassFermiFilter<T>,
    ComplementaryFilterAdaptor<LowPassFermiFilter<T>>>;

  template<typename T>
  using MultiLevelRecombinedFilterFamily = GenericMultiLevelFilterFamily<NullFilter<T>, Filter<T>>;

  template<typename T>
  using UnfilteredFilterFamily = GenericMultiLevelFilterFamily<Filter<T>, Filter<T>>;


  template<typename T>
  class ResidualFilterFamily : public FilterFamily<T> {

  public:
    ResidualFilterFamily(const FilterFamily<T> &underlying) {
      for (size_t i = 0; i < underlying.getMaxLevel() - 1; i++) {
        this->filters.push_back(underlying.complementFilters[i]->clone());
      }
      this->filters.push_back(std::make_shared<NullFilter<T>>());
    }

  };


  template<typename T>
  std::ostream &operator<<(std::ostream &s, const Filter<T> &f) {
    f.debugInfo(s);
    return s;
  }

  template<typename T>
  std::ostream &operator<<(std::ostream &s, const FilterFamily<T> &f) {
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
}

#endif //IC_FILTER_HPP
