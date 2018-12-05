#ifndef IC_FILTERFAMILY_HPP
#define IC_FILTERFAMILY_HPP

#include <memory>
#include <vector>
#include <stdexcept>
#include "src/simulation/filters/filter.hpp"

namespace filters {
  template<typename T>
  class ResidualFilterFamily;

  //! \class FilterFamily
  //! \brief Generic class to define filters on multiple levels of grids.
  template<typename T>
  class FilterFamily {
  protected:
    friend class ResidualFilterFamily<T>;

    std::vector<std::shared_ptr<filters::Filter<T>>> filters;//!< Filters for a given level
    std::vector<std::shared_ptr<filters::Filter<T>>> complementFilters; //!< Complementary filters for a given level
    std::vector<std::shared_ptr<filters::Filter<T>>> hpFilters;//!< High pass filters for a given level
    std::vector<std::shared_ptr<filters::Filter<T>>> lpFilters;//!< Low pass filters for a given level

    //! Default constructor
    FilterFamily() {

    }

  public:

    //! Constructor with known number of levels
    FilterFamily(size_t maxLevels) {
      for (size_t i = 0; i < maxLevels; ++i)
        filters.push_back(std::make_shared<Filter<T>>());
    }


    //! Returns the filter on the specified level
    virtual const Filter <T> &getFilterOnLevel(int level) const {
      return *(filters[level]);
    }

    //! Returns the high pass filter on the specified level
    virtual const Filter <T> &getHighPassFilterOnLevel(int level) const {
      return *(hpFilters[level]);
    }

    //! Returns the low pass filter on the specified level
    virtual const Filter <T> &getLowPassFilterOnLevel(int level) const {
      return *(lpFilters[level]);
    }

    //! Returns the number of levels of filter that have been defined.
    virtual size_t getMaxLevel() const {
      return filters.size();
    }

    //! Outputs debug information
    virtual void debugInfo(std::ostream &s) const {
      s << "FilterFamily(";
      for (size_t i = 0; i < filters.size(); ++i) {
        s << *filters[i];
        if (i + 1 < filters.size()) s << ", ";
      }
      s << ")";

    }

    //! Add a level to this filter family.
    virtual void addLevel(T /*k_cut*/) {
      throw std::runtime_error("Don't know how to add a level to this type of filter family");
    }
  };

  /*! \class GenericMultiLevelFilterFamily
    \brief Generic filter family that has a high-pass type and low-pass type filter
  */
  template<typename LowFilterType, typename HighFilterType>
  class GenericMultiLevelFilterFamily : public FilterFamily<typename LowFilterType::ReturnType> {
  protected:
    using T = typename LowFilterType::ReturnType;
    static_assert(std::is_same<typename LowFilterType::ReturnType, typename HighFilterType::ReturnType>::value,
                  "Filters must use same floating point type");
  public:

  //! Constructor - start with a single level of filters
    GenericMultiLevelFilterFamily() {
      this->filters.push_back(std::make_shared<Filter<T>>());
      this->lpFilters.push_back(std::make_shared<Filter<T>>());
      this->hpFilters.push_back(std::make_shared<Filter<T>>());
    }

    //! Adds relevant filters to the next level to be defined, based on what was used on the previous level
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
  using MultiLevelFilterFamily = GenericMultiLevelFilterFamily<filters::LowPassFermiFilter<T>,
      filters::ComplementaryCovarianceFilterAdaptor<filters::LowPassFermiFilter<T>>>;

  template<typename T>
  using MultiLevelDependentFilterFamily = GenericMultiLevelFilterFamily<filters::LowPassFermiFilter<T>,
      filters::ComplementaryFilterAdaptor<filters::LowPassFermiFilter<T>>>;

  template<typename T>
  using MultiLevelRecombinedFilterFamily = GenericMultiLevelFilterFamily<filters::NullFilter<T>, filters::Filter<T>>;

  template<typename T>
  using UnfilteredFilterFamily = GenericMultiLevelFilterFamily<filters::Filter<T>, filters::Filter<T>>;


  /*! \class ResidualFilterFamily
  \brief Copies another family of filters, but with the complement of each filter on each level.
  */
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


  //! Output debug information about a filter to a stream
  template<typename T>
  std::ostream &operator<<(std::ostream &s, const Filter <T> &f) {
    f.debugInfo(s);
    return s;
  }

    //! Output debug information about a family of filters to a stream
  template<typename T>
  std::ostream &operator<<(std::ostream &s, const FilterFamily<T> &f) {
    f.debugInfo(s);
    return s;
  }
}
#endif
