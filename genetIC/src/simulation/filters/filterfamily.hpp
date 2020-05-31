#ifndef IC_FILTERFAMILY_HPP
#define IC_FILTERFAMILY_HPP

#include <memory>
#include <vector>
#include <stdexcept>
#include "src/simulation/filters/filter.hpp"

namespace multilevelgrid {
  template<typename DataType, typename T>
  class MultiLevelGrid;
}

namespace filters {

  //! \class FilterFamilyBase
  //! \brief Generic class to define filters on multiple levels of grids.
  template<typename T>
  class FilterFamilyBase {
  protected:

    std::vector<std::shared_ptr<filters::Filter<T>>> filters;//!< Filters for a given level
    std::vector<std::shared_ptr<filters::Filter<T>>> complementFilters; //!< Complementary filters for a given level
    std::vector<std::shared_ptr<filters::Filter<T>>> hpFilters;//!< High pass filters for a given level
    std::vector<std::shared_ptr<filters::Filter<T>>> lpFilters;//!< Low pass filters for a given level

    //! Default constructor
    FilterFamilyBase() {

    }

  public:

    //! Constructor with known number of levels
    FilterFamilyBase(size_t maxLevels) {
      for (size_t i = 0; i < maxLevels; ++i)
        filters.push_back(std::make_shared<Filter<T>>());
    }


    //! Returns the filter on the specified level
    virtual const Filter <T> &getFilterForLevel(int level) const {
      return *(filters[level]);
    }

    //! Returns the high pass filter on the specified level
    virtual const Filter <T> &getHighPassFilterForLevel(int level) const {
      return *(hpFilters[level]);
    }

    //! Returns the low pass filter on the specified level
    virtual const Filter <T> &getLowPassFilterForLevel(int level) const {
      return *(lpFilters[level]);
    }

    //! Returns the number of levels of filter that have been defined.
    virtual size_t getMaxLevel() const {
      return filters.size();
    }

    //! Outputs debug information
    virtual void debugInfo(std::ostream &s) const {
      s << "FilterFamilyBase(";
      for (size_t i = 0; i < filters.size(); ++i) {
        s << *filters[i];
        if (i + 1 < filters.size()) s << ", ";
      }
      s << ")";

    }

    //! Add a level to this filter family.
    virtual void addLevel(T k_cut) = 0;
  };

  /*! \class FilterFamily
    \brief Implementation of default set of Fourier-space filters for combining information from fields on different levels
  */
  template<typename T>
  class FilterFamily : public FilterFamilyBase<T> {
  protected:
    using LowFilterType = filters::LowPassFermiFilter<T>; // this could be templatised to allow more flexibility
    using HighFilterType = filters::ComplementaryCovarianceFilterAdaptor<LowFilterType>;
  public:

    template<typename S>
    explicit FilterFamily(const multilevelgrid::MultiLevelGridBase<T, S> &fromContext) {
      if(fromContext.getLevelsAreCombined()) {
        // Output fields have been generated already...
        // Regard levels as independent. Strictly the filters should now be spatial, i.e. each grid
        // is used in its own domain of validity, but we assume we are only interested in the highest
        // resolution area (since in practice the only place the filters will now be used is in
        // checking modification values).
        auto nullfilter = std::make_shared<NullFilter<T>>();
        auto identityfilter = std::make_shared<Filter<T>>();
        for(size_t level=0; level < fromContext.getNumLevels() - 1; ++level) {
          this->filters.push_back(nullfilter);
          this->lpFilters.push_back(nullfilter);
          this->hpFilters.push_back(nullfilter);
        }
        this->filters.push_back(identityfilter);
        this->lpFilters.push_back(identityfilter);
        this->hpFilters.push_back(identityfilter);

      } else {
        this->filters.push_back(std::make_shared<Filter<T>>());
        this->lpFilters.push_back(std::make_shared<Filter<T>>());
        this->hpFilters.push_back(std::make_shared<Filter<T>>());

        for (size_t level = 0; level < fromContext.getNumLevels() - 1; ++level) {
          const grids::Grid<T> &grid0(fromContext.getGridForLevel(level));

          T k_pixel = ((T) grid0.size) * grid0.getFourierKmin();
          T k_cut = FRACTIONAL_K_SPLIT * k_pixel;
          this->addLevel(k_cut);
        }
      }
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


  //! Output debug information about a filter to a stream
  template<typename T>
  std::ostream &operator<<(std::ostream &s, const Filter <T> &f) {
    f.debugInfo(s);
    return s;
  }

  //! Output debug information about a family of filters to a stream
  template<typename T>
  std::ostream &operator<<(std::ostream &s, const FilterFamilyBase<T> &f) {
    f.debugInfo(s);
    return s;
  }
}
#endif
