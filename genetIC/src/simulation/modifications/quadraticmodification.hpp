#ifndef IC_QUADRATICMODIFICATION_HPP
#define IC_QUADRATICMODIFICATION_HPP

#include <numeric>
#include <functional>
#include <src/simulation/modifications/modification.hpp>

namespace modifications {


  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class QuadraticModification : public Modification<DataType, T> {
  protected:
    int initNumberSteps;
    T targetPrecision;    /*!< Precision at which the target should be achieved. 0.01 input is 1% required accuracy */

  public:
    QuadraticModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_) :
        Modification<DataType, T>(underlying_, cosmology_) {
      this->order = 2;
    };

    QuadraticModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_, int initNumberSteps_,
                          T targetPrecision_) :
        Modification<DataType, T>(underlying_, cosmology_), initNumberSteps(initNumberSteps_) {
      this->targetPrecision = targetPrecision_;
      this->order = 2;
    };

    int getInitNumberSteps() {
      return this->initNumberSteps;
    }

    T getTargetPrecision() {
      return this->targetPrecision;
    }


    T calculateCurrentValue(const fields::MultiLevelField<DataType> &field) override {

      auto pushedField = pushMultiLevelFieldThroughMatrix(field);
      pushedField->toFourier();
      T value = pushedField->innerProduct(field).real();
      return value;
    }

    std::shared_ptr<fields::ConstraintField<DataType>>
    pushMultiLevelFieldThroughMatrix(const fields::MultiLevelField<DataType> &field) {

      std::vector<std::shared_ptr<fields::Field<DataType, T>>> pushedfields;

      for (size_t level = 0; level < this->underlying.getNumLevels(); level ++){

        fields::Field<DataType, T> pushedField = this->pushOneLevelFieldThroughMatrix(field.getFieldForLevel(level), level);
        using tools::numerics::operator/=;
        pushedField.getDataVector() /= this->underlying.getWeightForLevel(level);
        pushedfields.push_back(std::make_shared<fields::Field<DataType, T>>(pushedField));
      }

      return std::make_shared<fields::ConstraintField<DataType>>(
          *dynamic_cast<multilevelcontext::MultiLevelContextInformation<DataType, T> *>(&(this->underlying)),
          pushedfields);
    }

  protected:
    virtual fields::Field<DataType, T>
    pushOneLevelFieldThroughMatrix(const fields::Field<DataType, T> &/* field */, size_t /* level */) = 0;
  };

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class FilteredVarianceModification : public QuadraticModification<DataType, T> {

  public:

    FilteredVarianceModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                                 const cosmology::CosmologicalParameters<T> &cosmology_, T filterscale_) :
        QuadraticModification<DataType, T>(underlying_, cosmology_), scale(filterscale_) {
      checkFilterScale(filterscale_);
      this->scale = filterscale_;
    }

    FilteredVarianceModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                                 const cosmology::CosmologicalParameters<T> &cosmology_, int initNumberSteps_,
                                 T targetPrecision_, T filterscale_) :
        QuadraticModification<DataType, T>(underlying_, cosmology_, initNumberSteps_, targetPrecision_){
      checkFilterScale(filterscale_);
      this->scale = filterscale_;
    }

    void checkFilterScale(T scale_) {

      size_t finest_level = this->underlying.getNumLevels() - 1;
      auto finest_grid = this->underlying.getGridForLevel(finest_level);
      T fine_pixelsize = finest_grid.cellSize;

      if (scale_ < fine_pixelsize) {
        throw std::runtime_error("Variance high-pass filtering scale is smaller than the smallest pixel size.");
      }

      T window_size = finest_grid.getFlaggedCellsPhysicalSize();
      if (scale_ > window_size) {
        std::cerr << "WARNING: High-pass filtering scale: " << scale_ << " h**-1 Mpc is greater than the rough window "
                  <<"size: " << window_size <<" h**-1 Mpc used for modifications." << std::endl;
        std::cerr << "This is prone to numerical errors when modifying the field."
                  << " Decrease filtering scale to avoid it." << std::endl;
      }

    }

    fields::Field<DataType, T> pushOneLevelFieldThroughMatrix(const fields::Field<DataType, T> &field, size_t level) override {

      fields::Field<DataType, T> pushedField = fields::Field<DataType, T>(field);

      pushedField.toReal();
      windowOperator(pushedField, level);

      pushedField.toFourier();
      filterOperator(pushedField);

      pushedField.toReal();
      varianceOperator(pushedField, level);

      pushedField.toFourier();
      filterOperator(pushedField);

      pushedField.toReal();
      windowOperator(pushedField, level);

      return pushedField;
    }

  private:
    T scale;

    void windowOperator(fields::Field<DataType, T> &field, size_t level) {

      assert(!field.isFourier()); // Windowing is done in real space

      std::vector<DataType> &fieldData = field.getDataVector();

#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < fieldData.size(); ++i) {
        // If cell is not a flagged cell, zero it
        if (!(std::binary_search(this->flaggedCells[level].begin(), this->flaggedCells[level].end(), i))) {
          fieldData[i] = 0;
        }
      }
    }

    void filterOperator(fields::Field<DataType, T> &field) {

      assert(field.isFourier());  //Filtering must be done in Fourier space

      T k_cut = 2 * M_PI / scale;

      auto highPassFermi = filters::ComplementaryFilterAdaptor<filters::LowPassFermiFilter<T>>(k_cut);

      field.applyFilter(highPassFermi);
    }

    void varianceOperator(fields::Field<DataType, T> &field, size_t level) {

      assert(!field.isFourier()); // Variance is calculated in real space

      windowOperator(field, level);

      std::vector<DataType> &fieldData = field.getDataVector();
      size_t regionSize = this->flaggedCells.size();

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
