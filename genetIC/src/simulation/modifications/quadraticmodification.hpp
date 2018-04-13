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

//      // Push the finest level and then down sample to other grids
//      size_t finest_level = this->underlying.getNumLevels() - 1;

      std::vector<std::shared_ptr<fields::Field<DataType, T>>> pushedfields;

      for (size_t level = 0; level < this->underlying.getNumLevels(); level ++){

        fields::Field<DataType, T> pushedField = this->pushOneLevelFieldThroughMatrix(field.getFieldForLevel(level));
        pushedfields.push_back(std::make_shared<fields::Field<DataType, T>>(pushedField));
      }

      return std::make_shared<fields::ConstraintField<DataType>>(
          *dynamic_cast<multilevelcontext::MultiLevelContextInformation<DataType, T> *>(&(this->underlying)),
          pushedfields);

//      using tools::numerics::operator/=;
//      auto highResPushedField = this->pushOneLevelFieldThroughMatrix(field.getFieldForLevel(level));
//      highResPushedField.toFourier();
//      highResPushedField.getDataVector() /= this->underlying.getWeightForLevel(level);
//      return this->underlying.generateMultilevelFromHighResField(std::move(highResPushedField));
    }

  protected:
    virtual fields::Field<DataType, T>
    pushOneLevelFieldThroughMatrix(const fields::Field<DataType, T> &/* field */) = 0;
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

//      size_t finest_level = this->underlying.getNumLevels() - 1;
//
//      if (finest_level == 0) {
//        T coarse_pixel = this->underlying.getGridForLevel(finest_level).cellSize;
//        if (scale_ > coarse_pixel) {
//          this->scale = scale_;
//        } else {
//          throw std::runtime_error("Variance filtering on scale smaller than pixel scale");
//        }
//      } else {
//        T coarse_pixel = this->underlying.getGridForLevel(finest_level - 1).cellSize;
//        T fine_pixel = this->underlying.getGridForLevel(finest_level).cellSize;
//        if (fine_pixel < scale_ && scale_ < (1.0 / 0.3) * coarse_pixel) {
//          this->scale = scale_;
//        } else {
//          throw std::runtime_error("Variance filtering scale must be kept far away from pixels scale");
//        }
//      }
    }

    fields::Field<DataType, T> pushOneLevelFieldThroughMatrix(const fields::Field<DataType, T> &field) override {

      fields::Field<DataType, T> pushedField = fields::Field<DataType, T>(field);

      pushedField.toReal();
      windowOperator(pushedField);

      pushedField.toFourier();
      filterOperator(pushedField);

      pushedField.toReal();
      varianceOperator(pushedField);

      pushedField.toFourier();
      filterOperator(pushedField);

      pushedField.toReal();
      windowOperator(pushedField);

      return pushedField;
    }

  private:
    T scale;

    void windowOperator(fields::Field<DataType, T> &field) {

      assert(!field.isFourier()); // Windowing is done in real space

      std::vector<DataType> &fieldData = field.getDataVector();

#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < fieldData.size(); ++i) {
        // If cell is not a flagged cell, zero it
        if (!(std::binary_search(this->flaggedCells.begin(), this->flaggedCells.end(), i))) {
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

    void varianceOperator(fields::Field<DataType, T> &field) {

      assert(!field.isFourier()); // Variance is calculated in real space

      windowOperator(field);

      std::vector<DataType> &fieldData = field.getDataVector();
      size_t regionSize = this->flaggedCells.size();

      // Calculate mean value in flagged region
      T sum = 0;
      for (size_t i = 0; i < regionSize; i++) {
        sum += fieldData[this->flaggedCells[i]];
      }

      for (size_t i = 0; i < regionSize; i++) {
        fieldData[this->flaggedCells[i]] -= sum / regionSize;
        fieldData[this->flaggedCells[i]] /= regionSize;
      }
    }

  };
}


#endif
