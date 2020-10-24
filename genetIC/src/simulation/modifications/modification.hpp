#ifndef IC_MODIFICATION_HPP
#define IC_MODIFICATION_HPP

#include "src/tools/data_types/complex.hpp"
#include "src/tools/logging.hpp"

namespace modifications {

  //! Abstract definition of a modification
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class Modification {
  private:
    T target;    /*!< Target to be achieved by the modification */

  protected:
    const multilevelgrid::MultiLevelGrid<DataType> &underlying; //!< Underlying multi-level context object.
    const cosmology::CosmologicalParameters<T> &cosmology; //!< Struct containing cosmological parameters.
    std::vector<std::vector<size_t>> flaggedCells; //!< Region targeted by the modification.
    unsigned int order; //!< Linear are first order, qudartic are second etc.
    particle::species forSpecies; //!< What type of output field are we applying this modification to?


  public:
    /*! \brief Constructor from underlying multi-level context and cosmological data.

          \param underlying_ - underlying multi-level context object.
          \param cosmology_ - struct containing cosmological parameters.
          \param forSpecies - specifies the nature of the field we will apply this modification to
      */
    Modification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                 const cosmology::CosmologicalParameters<T> &cosmology_) : underlying(underlying_),
                                                 cosmology(cosmology_),
                                                 flaggedCells(underlying_.getNumLevels()),
                                                 forSpecies(particle::species::unknown)  {
      for (size_t level = 0; level < this->underlying.getNumLevels(); level++) {
        auto grid = this->underlying.getGridForLevel(level);
        grid.getFlaggedCells(flaggedCells[level]);


        if (this->flaggedCells[level].size() == grid.size3 && level != 0) {
          logging::entry(logging::level::warning) << "WARNING: Region selected for modification is the entire zoom grid. This is likely "
                    << "because the cell selection extends beyond the zoom boundaries." << std::endl;
          std::cerr
            << "By design, modifications are meant to be defined inside a zoom region. Increase the size of your "
            << "zoom grid or decrease your selection to avoid nasty surprises." << std::endl;
        }
      }
    };

    //! Calculate modification value with a given field
    virtual T calculateCurrentValue(const fields::MultiLevelField<DataType> & /* field */) = 0;

    //! Returns the target of the modification, ie, what we want the function of the field to be constrained to be.
    T getTarget() {
      return target;
    }

    //! Sets the modification's target to the specified value.
    void setTarget(T target_) {
      target = target_;
    }

    //! Returns 1 for linear modifications, and 2 for quadratic modifications (and n for nth order modifications if these were created)
    unsigned int getOrder() {
      return this->order;
    }

  };


  //! \brief Exception to identify unknown modifications, if the user requests a modification that isn't recognised.
  class UnknownModificationException : public std::exception {
  public:

    explicit UnknownModificationException(const char *message_) : message(message_) {}

    explicit UnknownModificationException(const std::string &message_) : message(message_) {}

    virtual ~UnknownModificationException() throw() {}

    virtual const char *what() const throw() {
      return message.c_str();
    }

  protected:
    std::string message;
  };
};


#endif
