#ifndef IC_MODIFICATION_HPP
#define IC_MODIFICATION_HPP

#include <src/tools/data_types/complex.hpp>

namespace modifications {

  //! Abstract definition of a modification
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class Modification {
  private:
    T target;    /*!< Target to be achieved by the modification */

  protected:
    multilevelcontext::MultiLevelContextInformation<DataType> &underlying;
    const cosmology::CosmologicalParameters<T> &cosmology;
    std::vector<size_t> flaggedCells;    /*!< Region targeted by the modification */
    unsigned int order;                  /*!< Linear are first order, qudartic are second etc */


  public:
    Modification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                 const cosmology::CosmologicalParameters<T> &cosmology_) : underlying(underlying_),
                                                                           cosmology(cosmology_) {

      size_t finestlevel = this->underlying.getNumLevels() - 1;
      auto finestgrid = this->underlying.getGridForLevel(finestlevel);
      finestgrid.getFlaggedCells(flaggedCells);


      if(this->flaggedCells.size() == finestgrid.size3 && finestlevel != 0){
        std::cerr << "WARNING: Region selected for modification is the entire zoom grid. This is likely "
                  << "because the cell selection extends beyond the zoom boundaries. By design, "
                      "modifications are only defined inside the zoom region. Increase the size of your "
                      "zoom grid or decrease your selection to avoid nasty surprises. " << std::endl;
      }
    };

    //! Calculate modification value with a given field
    virtual T calculateCurrentValue(const fields::MultiLevelField<DataType> & /* field */) = 0;

    T getTarget() {
      return target;
    }

    void setTarget(T target_) {
      target = target_;
    }

    unsigned int getOrder() {
      return this->order;
    }

  };


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
