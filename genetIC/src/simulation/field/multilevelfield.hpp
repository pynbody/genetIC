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
    multilevelcontext::MultiLevelContextInformation<DataType> *multiLevelContext;
    std::shared_ptr<filters::FilterFamily<T>> pFilters;   // filters to be applied when used as a vector
    tools::Signaling::connection_t connection;
    bool isCovector;

    std::vector<Field<DataType, T>> fieldsOnLevels;

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
        T k_cut = ((T) grid0.size) * FRACTIONAL_K_SPLIT * 2. * M_PI / grid0.boxsize;
        this->pFilters->addLevel(k_cut);
      }

    }

    MultiLevelField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext) : multiLevelContext(&multiLevelContext) {
      setupConnection();
    }

    MultiLevelField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    std::vector<Field<DataType, T>> &&fieldsOnGrids) :
      multiLevelContext(&multiLevelContext), fieldsOnLevels(std::move(fieldsOnGrids)) {
      setupConnection();
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
      return fieldsOnLevels[i];
    }

    virtual const Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) const {
      return (const_cast<MultiLevelField<DataType> *>(this))->getFieldForGrid(grid);
    };

    virtual Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) {
      // TODO: problematically slow implementation
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        if (&grid == &multiLevelContext->getGridForLevel(i))
          return getFieldForLevel(i);
      }
      throw (std::runtime_error("Cannot find a field for the specified grid"));
    };

    virtual Field<DataType, T> &getFieldForLevel(size_t i) {
      return fieldsOnLevels[i];
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
      const filters::Filter<T> *pFiltThis;
      const grids::Grid<T> *pCurrentGrid;
      vector<DataType> *pFieldThis;


      auto newLevel = [&](size_t level) {
        pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
        pFiltThis = &(getFilterForLevel(level));
        pFieldThis = &(getFieldForLevel(level).getDataVector());
        return pFieldThis->size() > 0;
      };

      auto newCell = [&](size_t level, size_t cell, size_t cumu_i) {
        (*pFieldThis)[cell] /= ratio;
      };

      multiLevelContext->forEachCellOfEachLevel(newLevel, newCell);
    }

    void addScaled(const MultiLevelField &other, DataType scale) {
      const filters::Filter<T> *pFiltThis;
      const filters::Filter<T> *pFiltOther;
      const grids::Grid<T> *pCurrentGrid;
      vector<DataType> *pFieldThis;
      const vector<DataType> *pFieldOther;

      auto newLevel = [&](size_t level) {
        pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
        pFiltThis = &(getFilterForLevel(level));
        pFiltOther = &(other.getFilterForLevel(level));

        // adjust our field's fourier/real convention to match the other field
        if (other.getFieldForLevel(level).isFourier())
          getFieldForLevel(level).toFourier();
        else
          getFieldForLevel(level).toReal();

        // check we did that right!
        assert (getFieldForLevel(level).isFourier() == other.getFieldForLevel(level).isFourier());

        pFieldThis = &(getFieldForLevel(level).getDataVector());
        pFieldOther = &(other.getFieldForLevel(level).getDataVector());

        /* cerr << "addScaled " << level << ":";
        pFiltThis->debugInfo(cerr);
        cerr << " -> ";
        pFiltOther->debugInfo(cerr);
        cerr << endl;*/

        return this->hasFieldOnGrid(level) && other.hasFieldOnGrid(level);
      };

      auto newCell = [&](size_t level, size_t cell, size_t cumu_i) {
        T k_value = pCurrentGrid->getFourierCellAbsK(cell);
        T filt = (*pFiltOther)(k_value) / (*pFiltThis)(k_value);
        (*pFieldThis)[cell] += filt * (*pFieldOther)[cell] * scale;
      };

      multiLevelContext->forEachCellOfEachLevel(newLevel, newCell);
    }

    DataType innerProduct(const MultiLevelField<DataType> &other) const {

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

      bool covariance_weighted;

      if (other.isCovector)
        covariance_weighted = true;

      T weight;
      const filters::Filter<T> *pFiltOther;
      const grids::Grid<T> *pCurrentGrid;
      const Field<DataType> *pFieldThis;
      const std::vector<DataType> *pFieldDataThis;
      const std::vector<DataType> *pFieldDataOther;
      const std::vector<T> *pCov;

      auto newLevelCallback = [&](size_t level) {
        weight = multiLevelContext->getWeightForLevel(level);
        pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
        pCov = &(multiLevelContext->getCovariance(level));
        pFiltOther = &(other.getFilterForLevel(level));
        pFieldThis = &(this->getFieldForLevel(level));
        pFieldDataThis = &(this->getFieldForLevel(level).getDataVector());
        pFieldDataOther = &(other.getFieldForLevel(level).getDataVector());
        return pFieldDataOther->size() > 0 && pFieldDataThis->size() > 0;
      };

      auto getCellContribution = [&](size_t component, size_t i, size_t cumu_i) {
        T k_value = pCurrentGrid->getFourierCellAbsK(i);
        T inner_weight = weight * (*pFiltOther)(k_value);
        if (covariance_weighted) inner_weight *= ((*pCov)[i]) * weight;
        inner_weight *= tools::numerics::fourier::getFourierCellWeight(*pFieldThis, i);
        return inner_weight * std::conj((*pFieldDataThis)[i]) * (*pFieldDataOther)[i];
      };

      return multiLevelContext->accumulateOverEachCellOfEachLevel(newLevelCallback, getCellContribution);

    }

    void applyFilters() {
      const filters::Filter<T> *pFiltThis;
      const grids::Grid<T> *pCurrentGrid;
      std::vector<DataType> *pFieldThis;

      auto newLevel = [&](size_t level) {
        pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
        pFiltThis = &(this->getFilterForLevel(level));
        /*cerr << level << " ";
        pFiltThis->debugInfo(cerr);
        cerr <<  endl;*/
        pFieldThis = &(this->getFieldForLevel(level).getDataVector());
        return pFieldThis->size() > 0;
      };

      auto newCell = [&](size_t level, size_t cell, size_t cumu_i) {
        T k_value = pCurrentGrid->getFourierCellAbsK(cell);
        (*pFieldThis)[cell] *= (*pFiltThis)(k_value);
      };

      multiLevelContext->forEachCellOfEachLevel(newLevel, newCell);

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

      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {

        auto &field = getFieldForLevel(i).getDataVector();
        const auto &spectrum = multiLevelContext->getCovariance(i);
        T norm = T(field.size());

        for (size_t j = 0; j < field.size(); ++j) {
          if (spectrum[j] != 0)
            chi2 += pow(abs(field[j]), 2.0) / (spectrum[j] * norm);
        }
      }

      return chi2;

    }

  private:
    void applySpectrumOneGrid(std::vector<DataType> &field,
                              const std::vector<T> &spectrum,
                              const grids::Grid<T> &grid) {

#pragma omp parallel for
      for (size_t i = 0; i < grid.size3; i++) {
        field[i] *= sqrt(spectrum[i]);
      }
    }

    void multiplyByCovarianceOneGrid(std::vector<DataType> &field,
                                     const std::vector<T> &spectrum,
                                     const grids::Grid<T> &grid,
                                     T weight) {

#pragma omp parallel for
      for (size_t i = 0; i < grid.size3; i++) {
        field[i] *= weight * spectrum[i];
      }
    }

    void enforceSpectrumOneGrid(Field<DataType> &field,
                                const std::vector<T> &spectrum,
                                const grids::Grid<T> &grid) {
      T white_noise_norm = sqrt(T(grid.size3));

// #pragma omp parallel for
      for (size_t i = 0; i < grid.size3; i++) {
        int kx, ky, kz;
        std::tie(kx, ky, kz) = grid.getFourierCellCoordinate(i);
        complex<T> existingValue = tools::numerics::fourier::getFourierCoefficient(field, kx, ky, kz);
        T absExistingValue = abs(existingValue);
        T sqrt_spec = sqrt(spectrum[i]) * white_noise_norm;
        T a = sqrt_spec * existingValue.real() / absExistingValue;
        T b = sqrt_spec * existingValue.imag() / absExistingValue;

        tools::numerics::fourier::setFourierCoefficient(field, kx, ky, kz, a, b);
      }
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
        this->fieldsOnLevels.emplace_back(this->multiLevelContext->getGridForLevel(i), source.getFieldForLevel(i));
      }
      this->fieldsOnLevels.emplace_back(
        this->multiLevelContext->getGridForLevel(source.getContext().getNumLevels() - 1));

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

  public:
    OutputField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext) : MultiLevelField<DataType>(
      multiLevelContext) {
      outputState = PRE_SEPARATION;
    }

    void updateMultiLevelContext() override {
      assert(outputState == PRE_SEPARATION);
      this->template setupFilters<filters::MultiLevelFilterFamily<T>>();
      this->fieldsOnLevels.clear();
      this->multiLevelContext->forEachLevel([this](grids::Grid<T> &g) {
        this->fieldsOnLevels.emplace_back(g);
      });
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


  };


  template<typename DataType>
  class ConstraintField : public MultiLevelField<DataType> {

  protected:
    using T = typename MultiLevelField<DataType>::T;


  public:
    ConstraintField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    std::vector<Field<DataType, T>> &&fieldsOnGrids)
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
