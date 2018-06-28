#ifndef IC_FILTER_HPP
#define IC_FILTER_HPP

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

    Filter() {}

    Filter(T /*k_cut*/) {
      // provided for compatibility with filters that do take a wavenumber cut
    }

    virtual T operator()(T) const {
      return 1.0;
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

    NullFilter(T /*k_cut*/) {
      // provided for compatibility with filters that do take a wavenumber cut
    }

    T operator()(T /*x*/) const override {
      return 0.0;
    }

    std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<NullFilter<T>>();
    }

    virtual void debugInfo(std::ostream &s) const override {
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

    virtual T operator()(T x) const override {
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

    virtual T operator()(T x) const override {
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

    virtual std::shared_ptr<Filter<T>> clone() const override {
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

    virtual std::shared_ptr<Filter<T>> clone() const override {
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
          typeid(pOriginal) !=
          typeid(ComplementaryFilterAdaptor<UnderlyingType, T>)); // this is NOT a copy constructor!
      pUnderlying = pOriginal.clone();
    }


    T operator()(T k) const override {
      return 1. - (*pUnderlying)(k);
    }

    virtual std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<ComplementaryFilterAdaptor<UnderlyingType>>(*this);
    }

    virtual void debugInfo(std::ostream &s) const override {
      s << "ComplementaryFilterAdaptor(" << (*pUnderlying) << ")";
    }
  };
}

#endif //IC_FILTER_HPP
