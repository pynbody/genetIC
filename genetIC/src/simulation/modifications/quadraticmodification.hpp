#ifndef IC_QUADRATICMODIFICATION_HPP
#define IC_QUADRATICMODIFICATION_HPP

#include <numeric>
#include <functional>
#include "src/simulation/modifications/modification.hpp"
#include "src/tools/logging.hpp"

namespace modifications {


  /*! \class QuadraticModification
    \brief Defines quadratic modifications, ie, those defined with quadratic functions of the field.
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class QuadraticModification : public Modification<DataType, T> {
  protected:
    int initNumberSteps; //!< Initial number of steps required.
    T targetPrecision;    /*!< Precision at which the target should be achieved. 0.01 input is 1% required accuracy */

  public:
    //! Constructor that leaves initNumberSteps and targetPrecision unspecified
    QuadraticModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_) :
      Modification<DataType, T>(underlying_, cosmology_) {
      this->order = 2;
    };

    //! Constructor that specifies initNumberSteps and targetPrecision
    QuadraticModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_, int initNumberSteps_,
                          T targetPrecision_) :
      Modification<DataType, T>(underlying_, cosmology_), initNumberSteps(initNumberSteps_) {
      this->targetPrecision = targetPrecision_;
      this->order = 2;
    };

    //! Returns the initial number of steps
    int getInitNumberSteps() {
      return this->initNumberSteps;
    }

    //! Returns the target precision
    T getTargetPrecision() {
      return this->targetPrecision;
    }


    //! Calculates the current value of the quadratic modification evaluated on the field which we are trying to constrain
    T calculateCurrentValue(const fields::MultiLevelField<DataType> &field) override {

      auto pushedField = pushMultiLevelFieldThroughMatrix(field);
      pushedField->toFourier();
      T value = pushedField->innerProduct(field).real();
      return value;
    }

    //! Applies matrix operation to each level of the multi-level field supplied
    std::shared_ptr<fields::ConstraintField<DataType>>
    pushMultiLevelFieldThroughMatrix(const fields::MultiLevelField<DataType> &field) {

      std::vector<std::shared_ptr<fields::Field<DataType, T>>> pushedfields;

      for (size_t level = 0; level < this->underlying.getNumLevels(); level++) {
        pushedfields.push_back(std::make_shared<fields::Field<DataType, T>>(field.getFieldForLevel(level)));
      }

      auto result = std::make_shared<fields::ConstraintField<DataType>>(
        *dynamic_cast<const multilevelgrid::MultiLevelGrid<DataType, T> *>(&(this->underlying)),
        pushedfields, field.getTransferType(),
        false /* the fields still represent a vector, not a covector */ );

      result->applyTransferRatio(field.getTransferType(), particle::species::dm);
      // always want the variance of the dm overdensity field
      // The above applies C_DM^(1/2) if field.transferType is whitenoise. If field.transferType is already dm,
      // it's a null-op.

      for (size_t level = 0; level < this->underlying.getNumLevels(); level++) {
        this->pushOneLevelFieldThroughMatrix(result->getFieldForLevel(level), level);
      }

      // at this point we have Q | delta > (where matrices follow notation of Rey & Pontzen 2018, eq 28).
      // We now want to turn this into < delta | Q (and can assume Q to be Hermitian).
      result->convertToCovector();


      result->applyTransferRatio(field.getTransferType(), particle::species::dm);
      // the transfer ratio for the incoming vector; i.e. this is applied here because we will be forming
      // <result | field> to calculate the value of the quadratic functional

      return result;
    }

  protected:
    //! Applies matrix operation that defines the quadratic modification to the specified level of the multi-level field.
    virtual void
    pushOneLevelFieldThroughMatrix(fields::Field<DataType, T> &/* field */, size_t /* level */) = 0;

  };

  /*! \class FilteredVarianceModification
    \brief Modifies the variance of the field, applying a filter in Fourier space

  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class FilteredVarianceModification : public QuadraticModification<DataType, T> {

  public:

    //! Constructor to specify filter space
    FilteredVarianceModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                                 const cosmology::CosmologicalParameters<T> &cosmology_, T filterscale_) :
      QuadraticModification<DataType, T>(underlying_, cosmology_), scale(filterscale_) {
      checkFilterScale(filterscale_);
      this->scale = filterscale_;
    }

    //! Constructor to specify filter space in addition to initial number of steps and target precision.
    FilteredVarianceModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                                 const cosmology::CosmologicalParameters<T> &cosmology_, int initNumberSteps_,
                                 T targetPrecision_, T filterscale_) :
      QuadraticModification<DataType, T>(underlying_, cosmology_, initNumberSteps_, targetPrecision_) {
      checkFilterScale(filterscale_);
      this->scale = filterscale_;
    }

    //! Check to see if the user actually specified a filter scale, and if it is suitable if they did
    void checkFilterScale(T scale_) {
      if (scale_ < 0.0) {
        throw std::runtime_error("Trying to calculate filtered variance without initialising variance filtering scale."
                                 " Use filtering_scale command to do this.");
      }

      size_t finest_level = this->underlying.getNumLevels() - 1;
      auto finest_grid = this->underlying.getGridForLevel(finest_level);
      T fine_pixelsize = finest_grid.cellSize;

      if (scale_ < fine_pixelsize) {
        throw std::runtime_error("Variance high-pass filtering scale is smaller than the smallest pixel size.");
      }

      T window_size = finest_grid.getFlaggedCellsPhysicalSize();
      if (scale_ > window_size) {
        logging::entry(logging::level::warning) << "WARNING: High-pass filtering scale: " << scale_ << " h**-1 Mpc is greater than the rough window "
                  << "size: " << window_size << " h**-1 Mpc used for modifications." << std::endl;
        logging::entry() << "This is prone to numerical errors when modifying the field."
                  << " Decrease filtering scale to avoid it." << std::endl;
      }

    }

    void pushOneLevelFieldThroughMatrix(fields::Field<DataType, T> &field, size_t level) override {

      field.toReal();
      windowOperator(field, level);

      field.toFourier();
      filterOperator(field);

      field.toReal();
      varianceOperator(field, level);

      field.toFourier();
      filterOperator(field);

      field.toReal();
      windowOperator(field, level);
    }

  private:
    T scale; //!< Filter scale used by the modification

    //! Zeros out everything that has not been flagged.
    void windowOperator(fields::Field<DataType, T> &field, size_t level) {

      assert(!field.isFourier()); // Windowing is done in real space

      std::vector<DataType> &fieldData = field.getDataVector();

#pragma omp parallel for schedule(static) default(none) shared(fieldData, level)
      for (size_t i = 0; i < fieldData.size(); ++i) {
        // If cell is not a flagged cell, zero it
        if (!(std::binary_search(this->flaggedCells[level].begin(), this->flaggedCells[level].end(), i))) {
          fieldData[i] = 0;
        }
      }
    }

    //! Applies a high pass fermi filter to the given field
    void filterOperator(fields::Field<DataType, T> &field) {

      assert(field.isFourier());  //Filtering must be done in Fourier space

      T k_cut = 2 * M_PI / scale;

      auto highPassFermi = filters::ComplementaryFilterAdaptor<filters::LowPassFermiFilter<T>>(k_cut);

      field.applyFilter(highPassFermi);
    }

    //! Apply the variance operator to the field, such that field.varianceOperator(field) = variance(field)
    void varianceOperator(fields::Field<DataType, T> &field, size_t level) {

      assert(!field.isFourier()); // Variance is calculated in real space

      windowOperator(field, level);

      std::vector<DataType> &fieldData = field.getDataVector();
      size_t regionSize = this->flaggedCells[level].size();

      // Calculate mean value in flagged region
      T sum = 0;
      for (size_t i = 0; i < regionSize; i++) {
        sum += fieldData[this->flaggedCells[level][i]];
      }

      for (size_t i = 0; i < regionSize; i++) {
        fieldData[this->flaggedCells[level][i]] -= sum / regionSize;
        fieldData[this->flaggedCells[level][i]] /= regionSize;
      }
    }

  };
}


#endif
