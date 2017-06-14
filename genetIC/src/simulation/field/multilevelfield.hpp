#ifndef IC_MULTILEVELFIELD_HPP
#define IC_MULTILEVELFIELD_HPP

#include "src/simulation/multilevelcontext/multilevelcontext.hpp"
#include "src/simulation/filters/filterfamily.hpp"
#include "src/simulation/field/field.hpp"


namespace fields {

  /** Class to manage a field defined across multiple grids. */
  template<typename DataType>
  class MultiLevelField : public std::enable_shared_from_this<MultiLevelField<DataType>> {

  protected:
    using T = tools::datatypes::strip_complex<DataType>;
    using ComplexType = tools::datatypes::ensure_complex<DataType>;
    multilevelcontext::MultiLevelContextInformation<DataType> *multiLevelContext;
    std::shared_ptr<filters::FilterFamily<T>> pFilters;   // filters to be applied when used as a vector
    tools::Signaling::connection_t connection;
    bool isCovector;

    std::vector<std::shared_ptr<Field<DataType, T>>> fieldsOnLevels;

    void setupConnection() {
      connection = multiLevelContext->connect([this]() {
        this->updateMultiLevelContext();
      });
    }


  public:

    template<typename FilterType>
    void setupFilters() {
      const T FRACTIONAL_K_SPLIT = 0.3;

      size_t nLevels = this->multiLevelContext->getNumLevels();
      this->pFilters = make_shared<FilterType>();
      for (size_t level = 0; level < nLevels - 1; ++level) {
        const grids::Grid<T> &grid0(this->multiLevelContext->getGridForLevel(level));
        T k_cut = ((T) grid0.size) * FRACTIONAL_K_SPLIT * 2. * M_PI / grid0.thisGridSize;
        this->pFilters->addLevel(k_cut);
      }

    }

    MultiLevelField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext) : multiLevelContext(
      &multiLevelContext) {
      setupConnection();
      isCovector = false;
    }

    MultiLevelField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    const std::vector<std::shared_ptr<Field<DataType, T>>> &fieldsOnGrids) :
      multiLevelContext(&multiLevelContext), fieldsOnLevels(fieldsOnGrids) {
      setupConnection();
      isCovector = false;
    }

    virtual void updateMultiLevelContext() {

    }

    virtual multilevelcontext::MultiLevelContextInformation<DataType> &getContext() const {
      return const_cast<multilevelcontext::MultiLevelContextInformation<DataType> &>(*multiLevelContext);
    }

    const filters::FilterFamily<T> &getFilters() const {
      return *pFilters;
    }


    virtual const Field<DataType, T> &getFieldForLevel(size_t i) const {
      assert(i<fieldsOnLevels.size());
      return *(fieldsOnLevels[i]);
    }

    virtual const Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) const {
      return (const_cast<MultiLevelField<DataType> *>(this))->getFieldForGrid(grid);
    };

    virtual Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) {
      // TODO: problematically slow implementation
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        if(grid.pointsToGrid(&multiLevelContext->getGridForLevel(i)))
          return getFieldForLevel(i);
      }
      throw (std::runtime_error("Cannot find a field for the specified grid"));
    };

    virtual Field<DataType, T> &getFieldForLevel(size_t i) {
      return *(fieldsOnLevels[i]);
    }


    size_t getNumLevels() const {
      return multiLevelContext->getNumLevels();
    }

    bool hasFieldOnGrid(size_t i) const {
      return this->getFieldForLevel(i).getDataVector().size() > 0;
    }

    virtual const filters::Filter<T> &getFilterForLevel(size_t i) const {
      return pFilters->getFilterOnLevel(i);
    }

    virtual const filters::Filter<T> &getHighPassFilterForLevel(size_t i) const {
      return pFilters->getHighPassFilterOnLevel(i);
    }

    virtual const filters::Filter<T> &getLowPassFilterForLevel(size_t i) const {
      return pFilters->getLowPassFilterOnLevel(i);
    }

    void toReal() {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        getFieldForLevel(i).toReal();
    }

    void toFourier() {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        getFieldForLevel(i).toFourier();
    }

    bool isCompatible(const MultiLevelField<DataType> &other) const {
      return other.multiLevelContext == multiLevelContext;
    }

    bool isRealOnAllLevels() const {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        if (getFieldForLevel(i).isFourier()) return false;
      return true;
    }

    bool isFourierOnAllLevels() const {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        if (!getFieldForLevel(i).isFourier()) return false;
      return true;
    }


    void operator+=(const MultiLevelField<DataType> &other) {
      assert (isCompatible(other));
      addScaled(other, 1.0);
    }

    void operator/=(DataType ratio) {
      using namespace tools::numerics;

      for(size_t level=0; level<getNumLevels(); level++) {
        auto & data = getFieldForLevel(level).getDataVector();
        data/=ratio;
      }
    }

    //! Add a scaled multilevel field to the current one
    /*!
        The two fields can have different filters on each level,
        hence the need to compensate filter normalization.
     */
    void addScaled(const MultiLevelField &other, DataType scale) {
      assert(other.isFourierOnAllLevels());
      toFourier();

      for(size_t level=0; level<getNumLevels(); level++) {
        if(hasFieldOnGrid(level) && other.hasFieldOnGrid(level)) {
          Field<DataType> &fieldThis = getFieldForLevel(level);
          const Field<DataType> &fieldOther = other.getFieldForLevel(level);
          auto &filtOther = (other.getFilterForLevel(level));
          auto &filtThis = getFilterForLevel(level);
          T kMin = fieldThis.getGrid().getFourierKmin();
          fieldThis.forEachFourierCellInt([&fieldOther, &filtOther, &filtThis, kMin, scale]
                                              (ComplexType currentVal, int kx, int ky, int kz) {
            T k_value = kMin * sqrt(T(kx * kx) + T(ky * ky) + T(kz * kz));
            T filt = filtOther(k_value) / filtThis(k_value);
            return currentVal + scale*filt * fieldOther.getFourierCoefficient(kx, ky, kz);
          });
        }
      }

    }

		// TODO Most important function of the code. Should be documented step by step (Args are very hard to undertsand)
    ComplexType innerProduct(const MultiLevelField<DataType> &other) const {

      assert(isCompatible(other));
      if (!isCovector)
        throw (std::runtime_error(
          "The inner product can only be taken if one of the fields is regarded as a covector"));
      /*
       * To understand why this restriction is in place, see notes on 'covector approach to constraints'
       *
       * TODO: translate these notes into the paper or into inline documentation here
       * */

      assert(isFourierOnAllLevels() && other.isFourierOnAllLevels());
      // To take inner product with correct filters, we must have the fields in fourier space

      bool covariance_weighted = false;

      if (other.isCovector)
        covariance_weighted = true;

      T weight;
      const filters::Filter<T> *pFiltOther;
      const grids::Grid<T> *pCurrentGrid;
      const Field<DataType> *pFieldThis, *pFieldOther;
      const std::vector<DataType> *pFieldDataThis;
      const Field<DataType> *pCov;

      ComplexType result(0,0);

      for(size_t level=0; level<getNumLevels(); ++level) {
        weight = multiLevelContext->getWeightForLevel(level);
        pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
        pCov = &(multiLevelContext->getCovariance(level));
        pFiltOther = &(other.getFilterForLevel(level));
        pFieldThis = &(this->getFieldForLevel(level));
        pFieldDataThis = &(this->getFieldForLevel(level).getDataVector());
        pFieldOther = &(other.getFieldForLevel(level));
        T kMin = pCurrentGrid->getFourierKmin();
        if(pFieldOther!=nullptr && pFieldDataThis->size() > 0) {
          result+=pFieldThis->accumulateForEachFourierCell([&](tools::datatypes::ensure_complex<DataType> thisFieldVal,
                                                int kx, int ky, int kz) {
            auto otherFieldVal = pFieldOther->getFourierCoefficient(kx,ky,kz);
            T k_value = kMin*sqrt(T(kx)*T(kx)+T(ky)*T(ky)+T(kz)*T(kz));
            T inner_weight = weight * (*pFiltOther)(k_value);
            if (covariance_weighted) inner_weight *= (pCov->getFourierCoefficient(kx,ky,kz).real()) * weight;
            return inner_weight * std::real(std::conj(thisFieldVal)*otherFieldVal);
          });
        }
      }
      return result;
    }

    T euclidianInnerProduct(const MultiLevelField<DataType> &other){
      // TODO implement this with std library
    }

    void applyFilters() {
      for(size_t level=0; level<getNumLevels(); ++level) {
        if(hasFieldOnGrid(level)) {
          getFieldForLevel(level).applyFilter(getFilterForLevel(level));
        }
      }

      pFilters = make_shared<filters::FilterFamily<T>>(multiLevelContext->getNumLevels());
    }


    void convertToVector() {
      assert(isCovector);
      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);

        multiplyByCovarianceOneGrid(getFieldForLevel(i),
                                    multiLevelContext->getCovariance(i),
                                    grid,
                                    multiLevelContext->getWeightForLevel(i));

      }
      isCovector = false;
    }

    void applyPowerSpectrum() {
      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);

        applySpectrumOneGrid(getFieldForLevel(i),
                             multiLevelContext->getCovariance(i),
                             grid);

      }
    }

    void enforceExactPowerSpectrum() {
      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);

        enforceSpectrumOneGrid(getFieldForLevel(i),
                               multiLevelContext->getCovariance(i),
                               grid);

      }
    }

    T getChi2() {

      T chi2 = 0;


      return chi2;

    }

  private:
    void applySpectrumOneGrid(Field<DataType> &field,
                              const Field<DataType> &spectrum,
                              const grids::Grid<T> &grid) {

      field.forEachFourierCellInt([&grid, &field, &spectrum]
                                   (complex<T> existingValue, int kx, int ky, int kz) {
        T sqrt_spec = sqrt(spectrum.getFourierCoefficient(kx,ky,kz).real());
        return existingValue*sqrt_spec;
      });
    }

    void multiplyByCovarianceOneGrid(Field<DataType> &field,
                                     const Field<DataType> &spectrum,
                                     const grids::Grid<T> &grid,
                                     T weight) {

      field.forEachFourierCellInt([weight, &grid, &field, &spectrum]
                                   (complex<T> existingValue, int kx, int ky, int kz) {
        T spec = spectrum.getFourierCoefficient(kx,ky,kz).real() * weight;
        return existingValue*spec;
      });
    }

    void enforceSpectrumOneGrid(Field<DataType> &field,
                                const Field<DataType> &spectrum,
                                const grids::Grid<T> &grid) {
      T white_noise_norm = sqrt(T(grid.size3));

      field.forEachFourierCell([white_noise_norm, &grid, &field, &spectrum]
                                   (complex<T> existingValue, int kx, int ky, int kz) {
        T absExistingValue = abs(existingValue);
        T sqrt_spec = sqrt(spectrum.getFourierCoefficient(kx,ky,kz).real()) * white_noise_norm;
        T a = sqrt_spec * existingValue.real() / absExistingValue;
        T b = sqrt_spec * existingValue.imag() / absExistingValue;
        return complex<T>(a,b);
      });

    }


  };


  template<typename DataType>
  class ResidualField : public MultiLevelField<DataType> {
  protected:
    using T = typename MultiLevelField<DataType>::T;

  public:
    ResidualField(const MultiLevelField<DataType> &source)
      : MultiLevelField<DataType>(source.getContext()) {
      for (size_t i = 0; i < source.getContext().getNumLevels() - 1; ++i) {
        this->fieldsOnLevels.emplace_back(
            std::make_shared<Field<DataType, T>>(this->multiLevelContext->getGridForLevel(i), source.getFieldForLevel(i))
        );
      }
      this->fieldsOnLevels.emplace_back(
        std::make_shared<Field<DataType, T>>(this->multiLevelContext->getGridForLevel(source.getContext().getNumLevels() - 1))
        );

      this->pFilters = std::make_shared<filters::ResidualFilterFamily<T>>(source.getFilters());

    }


    virtual void updateMultiLevelContext() override {
      throw (std::runtime_error("Update to grid structure took place while a ResidualField was in scope"));
    }

  };


  template<typename DataType>
  class OutputField : public MultiLevelField<DataType> {

  protected:
    using T = typename MultiLevelField<DataType>::T;
    typedef enum {
      PRE_SEPARATION,
      SEPARATED,
      RECOMBINED
    } t_output_state;

    t_output_state outputState;
    bool fieldsOnLevelsPopulated;

    void populateFieldsOnLevels() {
      this->multiLevelContext->forEachLevel([this](grids::Grid<T> &g) {
        this->fieldsOnLevels.emplace_back(std::make_shared<Field<DataType, T>>(g));
      });
      fieldsOnLevelsPopulated=true;
    }

    void populateFieldsOnLevelsIfRequired() {
      if(!fieldsOnLevelsPopulated)
        populateFieldsOnLevels();
    }

  public:
    OutputField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext)
      : MultiLevelField<DataType>(
      multiLevelContext) {
      outputState = PRE_SEPARATION;
      fieldsOnLevelsPopulated=false;
    }

    void updateMultiLevelContext() override {
      assert(outputState == PRE_SEPARATION);
      this->template setupFilters<filters::MultiLevelFilterFamily<T>>();
      this->fieldsOnLevels.clear();
      fieldsOnLevelsPopulated=false;
    }

    auto getHighKResiduals() {
      assert(outputState == PRE_SEPARATION);
      outputState = SEPARATED;
      return ResidualField<DataType>(*this);
    }

    void recombineHighKResiduals(const ResidualField<DataType> &residuals) {
      assert(outputState == SEPARATED);
      outputState = RECOMBINED;
      (*this) += residuals;
      this->template setupFilters<filters::MultiLevelRecombinedFilterFamily<T>>();
    }

    void setStateRecombined() {
      outputState = RECOMBINED;
      this->template setupFilters<filters::MultiLevelRecombinedFilterFamily<T>>();
    }

    Field<DataType, T> &getFieldForLevel(size_t i) override {
      populateFieldsOnLevelsIfRequired();
      return *(this->fieldsOnLevels[i]);
    }


  };


  template<typename DataType>
  class ConstraintField : public MultiLevelField<DataType> {

  protected:
    using T = typename MultiLevelField<DataType>::T;


  public:
    ConstraintField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    const std::vector<std::shared_ptr<Field<DataType, T>>> &fieldsOnGrids)
      : MultiLevelField<DataType>(multiLevelContext, std::move(fieldsOnGrids)) {
      this->isCovector = true;
      updateMultiLevelContext();
    }


    virtual void updateMultiLevelContext() override {
      this->template setupFilters<filters::MultiLevelDependentFilterFamily<T>>();
    }


  };
}

#endif //IC_MULTILEVELFIELD_HPP
