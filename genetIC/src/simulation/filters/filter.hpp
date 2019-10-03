#ifndef IC_FILTER_HPP
#define IC_FILTER_HPP

#include <stdexcept>
/*!
    \namespace filters
    \brief Define the filters used to separate low and high k modes in the box

 */
namespace filters {

  template<typename T>
  class ValueAsFilter;

  template<typename T>
  class DivisionFilter;

  template<typename T>
  class ProductFilter;

  /*! \class Filter
    \brief Class which is used to suppress the field at some values of k, and not at others.

    Essential for implementing the code's handling of zoom regions, as we lack access the high resolution
    data in most of the simulation, and have to find a way of cutting off fields without introducing too hard
    a cutoff.

  */
  template<typename T>
  class Filter : public std::enable_shared_from_this<Filter<T>> {
  public:
    typedef T ReturnType;

    //! Default constructor
    Filter() {}

    //! \brief Returns the filter at the specified k
    /*!
    \param - k - value of k at which to evaluate the filter
    */
    virtual T operator()(T) const {
      return 1.0;
    }

    //! Returns a copy of the filter
    virtual std::shared_ptr<Filter<T>> clone() const {
      return std::make_shared<Filter<T>>();
    }

    //! Returns a division with this filter as a numerator, the other as the denominator
    DivisionFilter<T> operator/(const Filter<T> &other) const {
      return DivisionFilter<T>(*this, other);
    }

    //! Returns a division with this filter as a numerator and a constant value on the denominator
    DivisionFilter<T> operator/(T value) const {
      return DivisionFilter<T>(*this, ValueAsFilter<T>(value));
    }

    //! Returns a product filter with this filter as one of the product members
    ProductFilter<T> operator*(const Filter<T> &other) const {
      return ProductFilter<T>(*this, other);
    }

    //! Returns a product between this filter and a constant value
    ProductFilter<T> operator*(T value) const {
      return ProductFilter<T>(*this, ValueAsFilter<T>(value));
    }

    //! Output debug information (debug use only)
    virtual void debugInfo(std::ostream &s) const {
      s << "Filter";
    }
  };

  //! \class NullFilter
  //! \brief Filters out everything (0 everywhere)
  template<typename T>
  class NullFilter : public Filter<T> {
  public:

    //! Default constructor
    NullFilter() {
    }


    //! Returns 0. Everything is filtered out by the null filter.
    T operator()(T /*x*/) const override {
      return 0.0;
    }

    //! Returns a copy of this filter.
    std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<NullFilter<T>>();
    }

    //! Output debug information to the specified stream.
    virtual void debugInfo(std::ostream &s) const override {
      s << "NullFilter";
    }

  };

  //! \class ValueAsFilter
  //! \brief Returns a single number regardless of the k value
  template<typename T>
  class ValueAsFilter : public Filter<T> {
  protected:
    T value;
  public:
    explicit ValueAsFilter(T value) : value(value) {
    }

    //! Returns the constant value
    T operator()(T /*x*/) const override {
      return value;
    }

    //! Returns a copy of this filter.
    std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<ValueAsFilter<T>>(value);
    }

    //! Output debug information to the specified stream.
    virtual void debugInfo(std::ostream &s) const override {
      s << value;
    }

  };


  /*! \class DivisionFilter
    \brief Filter obtained by the ratio of two other filters
  */
  template<typename T>
  class DivisionFilter : public Filter<T> {
  protected:
    std::shared_ptr<Filter<T>> pFirst; //!< Pointer to the first filter (numerator)
    std::shared_ptr<Filter<T>> pSecond; //!< Pointer to the second filter (denominator)
  public:
    //! \brief Constructor
    /*!
    \param first - Numerator filter
    \param second - Denominator filter
    */
    DivisionFilter(const Filter<T> &first, const Filter<T> &second) : pFirst(first.clone()), pSecond(second.clone()) {

    }

    //! Return the ratio of the two filters evaluated at the specified point
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

  /*! \class ProductFilter
    \brief A filter formed from the product of two other filters
  */
  template<typename T>
  class ProductFilter : public DivisionFilter<T> {
  protected:
    using DivisionFilter<T>::pFirst;
    using DivisionFilter<T>::pSecond;

  public:

    //! Constructor - takes the two filters of the product as arguments
    ProductFilter(const Filter<T> &first, const Filter<T> &second) : DivisionFilter<T>(first, second) {

    }

    //! Clone this filter and its sub-filters.
    virtual std::shared_ptr<Filter<T>> clone() const override {
      // simplification/optimization:
      if (typeid(*pFirst) == typeid(Filter<T>))
        return pSecond->clone();
      else if (typeid(*pSecond) == typeid(Filter<T>))
        return pFirst->clone();
      else
        return std::make_shared<ProductFilter<T>>(*pFirst, *pSecond);
    }

    //! Return the product of the two filters
    virtual T operator()(T x) const override {
      return (*pFirst)(x) * (*pSecond)(x);
    }

    //! Output debug information (debug use only)
    virtual void debugInfo(std::ostream &s) const override {
      s << (*pFirst) << "*" << (*pSecond);
    }
  };


  /*!
  \class LowPassFermiFilter
  \brief Low-pass filter constructed from the Fermi-dirac distribution
  'Temperature' is hard-coded at 1/10 of the cutoff value.
  */
  template<typename T>
  class LowPassFermiFilter : public Filter<T> {
  private:
    T kcut; //!< Cutoff value
    T temperature; //!< 'Temperature' of the fermi dirac distribution

  public:
    //! Constructor, with cutoff as argument
    LowPassFermiFilter(T kcut_in) : kcut(kcut_in) {
      // The following arbitrary choice for the temperature (i.e.
      // sharpness of the cut-off) has been made by numerical experiments
      // trading off the various pros and cons of a sharp cut discussed
      // in the notes.
      temperature = kcut * FRACTIONAL_FILTER_TEMPERATURE;
    };


    //! Return a copy of this filter
    virtual std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<LowPassFermiFilter<T>>(*this);
    }

    //! Fermi dirac distribution as the filter's temperature and cutoff in k-space
    T operator()(T k) const override {
      return 1. / (1. + exp((k - kcut) / temperature));
    }

    //! Debug information (only used for debugging purposes)
    virtual void debugInfo(std::ostream &s) const override {
      s << "LowPassFermiFilter(" << kcut << ")";
    }
  };


  /*! \class ComplementaryCovarianceFilterAdaptor
  \brief Complementary filter when applied to covariance
  Returns sqrt{1 - f(k)^2} for filter f
  */
  template<typename UnderlyingType, typename T=typename UnderlyingType::ReturnType>
  class ComplementaryCovarianceFilterAdaptor : public Filter<T> {
  private:

    std::shared_ptr<Filter<T>> pUnderlying;
  public:

    //! Constructor - arguments depend on the filter we are complementing
    template<typename... Args>
    ComplementaryCovarianceFilterAdaptor(Args &&... args) {
      pUnderlying = std::make_shared<UnderlyingType>(std::forward<Args>(args)...);
    };


    //! Returns \sqrt{1 - f(k)^2} for filter f
    T operator()(T k) const override {
      T denFilter = (*pUnderlying)(k);
      return sqrt(1. - denFilter * denFilter);
    }

    //! Create a copy of the filter
    virtual std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<ComplementaryCovarianceFilterAdaptor<UnderlyingType>>(*this);
    }

    //! Debug information (only used for debugging purposes)
    virtual void debugInfo(std::ostream &s) const override {
      s << "ComplementaryCovarianceFilterAdaptor(" << (*pUnderlying) << ")";
    }
  };

  /*! \class ComplementaryFilterAdaptor
      \brief Naive complementary filter

       This filter is constructed from another filter f, and returns 1 - f(k).
  */
  template<typename UnderlyingType, typename T=typename UnderlyingType::ReturnType>
  class ComplementaryFilterAdaptor : public Filter<T> {
  private:
    std::shared_ptr<Filter<T>> pUnderlying; //!< Pointer to the underlying filter.
  public:

    //! Constructor taking arbitrary arguments, which are used to construct the underlying filter
    template<typename... Args>
    ComplementaryFilterAdaptor(Args &&... args) {
      pUnderlying = std::make_shared<UnderlyingType>(std::forward<Args>(args)...);
    };

    //! Construct the complement to an existing filter (i.e. 1-exiting_filter)
    ComplementaryFilterAdaptor(const Filter<T> &pOriginal) {
      assert(
        typeid(pOriginal) !=
        typeid(ComplementaryFilterAdaptor<UnderlyingType, T>)); // this is NOT a copy constructor!
      pUnderlying = pOriginal.clone();
    }


    //!Returns 1 - f(k) for filter f
    T operator()(T k) const override {
      return 1. - (*pUnderlying)(k);
    }

    //! Create a copy of this filter
    virtual std::shared_ptr<Filter<T>> clone() const override {
      return std::make_shared<ComplementaryFilterAdaptor<UnderlyingType>>(*this);
    }

    //! Debug information (only used for debugging purposes)
    virtual void debugInfo(std::ostream &s) const override {
      s << "ComplementaryFilterAdaptor(" << (*pUnderlying) << ")";
    }
  };
}

#endif //IC_FILTER_HPP
