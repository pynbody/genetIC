#ifndef IC_MULTILEVELFIELD_HPP
#define IC_MULTILEVELFIELD_HPP

#include "multilevelcontext.hpp"
#include "filter.hpp"

template<typename T>
class strip_complex_s {
public:
  using type = T;
};

template<typename T>
class strip_complex_s<std::complex<T>> {
public:
  using type = T;
};


template <typename T>
using strip_complex = typename strip_complex_s<T>::type;




template<typename DataType>
class MultiLevelField {

protected:
  using T = strip_complex<DataType>;
  MultiLevelContextInformation<T> * multiLevelContext;
  std::shared_ptr<FilterFamily<T>> pFilters;   // filters to be applied when used as a vector
  Signaling::connection_t connection;
  bool isCovector;

  template<typename FilterType>
  void setupFilters() {
    size_t nLevels = this->multiLevelContext->getNumLevels();
    if (nLevels == 2) {
      const T FRACTIONAL_K_SPLIT = 0.3;
      const Grid<T> & grid0(this->multiLevelContext->getGridForLevel(0));
      this->pFilters = make_shared<FilterType>(
              ((T) grid0.size) * FRACTIONAL_K_SPLIT * 2. * M_PI / grid0.boxsize);
    } else if (nLevels == 1) {
      this->pFilters = make_shared<FilterFamily<T>>(1);
    } else {
      this->pFilters = nullptr;
      // throw (runtime_error("Don't know how to make appropriate filters for multi-level field"));
    }
  }


public:

  MultiLevelField(MultiLevelContextInformation<T> & multiLevelContext) : multiLevelContext(&multiLevelContext) {
    connection = multiLevelContext.connect([this]() {
      this->updateMultiLevelContext();
    });
  }

  virtual void updateMultiLevelContext() {

  }

  virtual MultiLevelContextInformation<T> & getContext() const {
    return const_cast<MultiLevelContextInformation<T> &>(*multiLevelContext);
  }

  const FilterFamily<T> & getFilters() const {
    return *pFilters;
  }

  virtual const std::vector<DataType> & getFieldOnGrid(size_t i) const = 0;

  std::vector<DataType> & getFieldOnGrid(size_t i) {
    auto & vec = const_cast<const MultiLevelField *>(this)->getFieldOnGrid(i);
    return const_cast<std::vector<DataType> &>(vec);
  }

  bool hasFieldOnGrid(size_t i) const {
    return this->getFieldOnGrid(i).size()>0;
  }

  virtual const Filter<T> & getFilterForGrid(size_t i) const {
    return pFilters->getFilterOnLevel(i);
  }

  bool isCompatible(const MultiLevelField<DataType> &other) const {
    return other.multiLevelContext == multiLevelContext;
  }


  void operator+=(const MultiLevelField<DataType> &other) {
    assert (isCompatible(other));
    addScaled(other, 1.0);
  }

  void operator/=(DataType ratio) {
    const Filter<T> *pFiltThis;
    const Grid<T> *pCurrentGrid;
    vector<DataType> *pFieldThis;


    auto newLevel = [&](size_t level) {
      pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
      pFiltThis = &(getFilterForGrid(level));
      pFieldThis = &(getFieldOnGrid(level));
      return pFieldThis->size()>0;
    };

    auto newCell = [&](size_t level, size_t cell , size_t cumu_i, vector<DataType> & unknown_field) {
      (*pFieldThis)[cell]/=ratio;
    };

    multiLevelContext->forEachCellOfEachLevel(newLevel, newCell);
  }

  void addScaled(const MultiLevelField &other, DataType scale) {
    const Filter<T> *pFiltThis;
    const Filter<T> *pFiltOther;
    const Grid<T> *pCurrentGrid;
    vector<DataType> *pFieldThis;
    const vector<DataType> *pFieldOther;

    auto newLevel = [&](size_t level) {
      pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
      pFiltThis = &(getFilterForGrid(level));
      pFiltOther = &(other.getFilterForGrid(level));
      pFieldThis = &(getFieldOnGrid(level));
      pFieldOther = &(other.getFieldOnGrid(level));
      return this->hasFieldOnGrid(level) && other.hasFieldOnGrid(level);
    };

    auto newCell = [&](size_t level, size_t cell , size_t cumu_i, vector<DataType> & unknown_field) {
      T k_value = pCurrentGrid->getFourierCellAbsK(cell);
      T filt = (*pFiltOther)(k_value)/(*pFiltThis)(k_value);
      (*pFieldThis)[cell] += filt * (*pFieldOther)[cell] * scale;
    };

    multiLevelContext->forEachCellOfEachLevel(newLevel, newCell);
  }

  DataType innerProduct(const MultiLevelField<DataType> &other) const {

    assert(isCompatible(other));
    if(!isCovector)
      throw(std::runtime_error("The inner product can only be taken if one of the fields is regarded as a covector"));
    /*
     * To understand why this restriction is in place, see notes on 'covector approach to constraints'
     *
     * TODO: translate these notes into the paper or into inline documentation here
     * */

    bool covariance_weighted;

    if(other.isCovector)
      covariance_weighted = true;

    T weight;
    const Filter<T> *pFiltOther;
    const Grid<T> *pCurrentGrid;
    const std::vector<DataType> *pFieldThis;
    const std::vector<DataType> *pFieldOther;
    const std::vector<T> *pCov;

    auto newLevelCallback = [&](size_t level) {
      weight = multiLevelContext->getWeightForLevel(level);
      pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
      pCov = &(multiLevelContext->getCovariance(level));
      pFiltOther = &(other.getFilterForGrid(level));
      pFieldThis = &(this->getFieldOnGrid(level));
      pFieldOther = &(other.getFieldOnGrid(level));
      return pFieldOther->size()>0 && pFieldThis->size()>0;
    };

    auto getCellContribution = [&](size_t component, size_t i, size_t cumu_i) {
      T k_value = pCurrentGrid->getFourierCellAbsK(i);
      T inner_weight = weight * (*pFiltOther)(k_value);
      if(covariance_weighted)  inner_weight*=((*pCov)[i])*weight;
      return inner_weight*conj((*pFieldThis)[i]) * (*pFieldOther)[i];
    };

    return multiLevelContext->accumulateOverEachCellOfEachLevel(newLevelCallback, getCellContribution);

  }

  void applyFilters() {
    const Filter<T> *pFiltThis;
    const Grid<T> *pCurrentGrid;
    std::vector<DataType> *pFieldThis;


    auto newLevel = [&](size_t level) {
      pCurrentGrid = &(multiLevelContext->getGridForLevel(level));
      pFiltThis = &(this->getFilterForGrid(level));
      pFieldThis = &(this->getFieldOnGrid(level));
      return pFieldThis->size()>0;
    };

    auto newCell = [&](size_t level, size_t cell , size_t cumu_i, std::vector<DataType> & unknown_field) {
      T k_value = pCurrentGrid->getFourierCellAbsK(cell);
      (*pFieldThis)[cell] *= (*pFiltThis)(k_value);
    };

    multiLevelContext->forEachCellOfEachLevel(newLevel, newCell);

    pFilters = make_shared<FilterFamily<T>>(multiLevelContext->getNumLevels());
  }



  void convertToVector() {
    assert(isCovector);
    for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
      auto & grid = multiLevelContext->getGridForLevel(i);

      multiplyByCovarianceOneGrid(getFieldOnGrid(i),
                                  multiLevelContext->getCovariance(i),
                                  grid,
                                  multiLevelContext->getWeightForLevel(i));

    }
    isCovector = false;
  }

  void applyPowerSpectrum() {
    for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
      auto & grid = multiLevelContext->getGridForLevel(i);

      applySpectrumOneGrid(getFieldOnGrid(i),
                           multiLevelContext->getCovariance(i),
                           grid);

    }
  }

  void enforceExactPowerSpectrum() {
    for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
      auto & grid = multiLevelContext->getGridForLevel(i);

      enforceSpectrumOneGrid(getFieldOnGrid(i),
                           multiLevelContext->getCovariance(i),
                           grid);

    }
  }

private:
  void applySpectrumOneGrid(std::vector<std::complex<T>> &field,
                            const std::vector<T> &spectrum,
                            const Grid<T> &grid) {

#pragma omp parallel for
    for (size_t i = 0; i < grid.size3; i++) {
      field[i] *= sqrt(spectrum[i]);
    }
  }

  void multiplyByCovarianceOneGrid(std::vector<std::complex<T>> &field,
                            const std::vector<T> &spectrum,
                            const Grid<T> &grid,
                            T weight)
  {

#pragma omp parallel for
    for (size_t i = 0; i < grid.size3; i++) {
      field[i] *= weight*spectrum[i];
    }
  }

  void enforceSpectrumOneGrid(std::vector<std::complex<T>> &field,
                            const std::vector<T> &spectrum,
                            const Grid<T> &grid) {
    T white_noise_norm = sqrt(T(grid.size3));

#pragma omp parallel for
    for (size_t i = 0; i < grid.size3; i++) {
      T existing_norm = abs(field[i]);
      field[i] *= sqrt(spectrum[i]) * white_noise_norm / existing_norm;
    }
  }


};



template<typename DataType>
class ResidualField : public MultiLevelField<DataType> {
protected:
  using T = typename MultiLevelField<DataType>::T;
  std::vector<std::vector<DataType>> fieldsOnGrids;

public:
  ResidualField(const MultiLevelField<DataType> & source)
    : MultiLevelField<DataType>(source.getContext())
  {
    for(size_t i=0; i<source.getContext().getNumLevels()-1; ++i) {
      fieldsOnGrids.push_back(source.getFieldOnGrid(i));
    }
    fieldsOnGrids.emplace_back();

    this->pFilters = std::make_shared<ResidualFilterFamily<T>>(source.getFilters());

  }


  virtual void updateMultiLevelContext() override {
    throw(std::runtime_error("Update to grid structure took place while a ResidualField was in scope"));
  }

  virtual const std::vector<DataType> &getFieldOnGrid(size_t i) const override {
    return fieldsOnGrids[i];
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
  OutputField(MultiLevelContextInformation<T> & multiLevelContext) : MultiLevelField<DataType>(multiLevelContext) {
    outputState = PRE_SEPARATION;
  }

  virtual void updateMultiLevelContext() {
    assert(outputState==PRE_SEPARATION);
    this->template setupFilters<TwoLevelFilterFamily<T>>();
  }

  auto getHighKResiduals() {
    assert(outputState==PRE_SEPARATION);
    outputState = SEPARATED;
    return ResidualField<DataType>(*this);
  }

  void recombineHighKResiduals(const ResidualField<DataType> &residuals) {
    assert(outputState==SEPARATED);
    outputState = RECOMBINED;
    (*this)+=residuals;
    this->template setupFilters<TwoLevelRecombinedFilterFamily<T>>();
  }


  virtual const std::vector<DataType> &getFieldOnGrid(size_t i) const override {
    return this->multiLevelContext->getGridForLevel(i).getFieldFourier();
  }
};



template<typename DataType>
class ConstraintField : public MultiLevelField<DataType> {

protected:
  using T = typename MultiLevelField<DataType>::T;
  std::vector<std::vector<DataType>> fieldsOnGrids;

public:
  ConstraintField(MultiLevelContextInformation<T> & multiLevelContext, std::vector<std::vector<DataType>> && fieldsOnGrids)
  : MultiLevelField<DataType>(multiLevelContext), fieldsOnGrids(std::move(fieldsOnGrids))
  {
    this->isCovector = true;
    updateMultiLevelContext();
  }


  virtual void updateMultiLevelContext() override {
    this->template setupFilters<TwoLevelDependentFilterFamily<T>>();
  }

  virtual const std::vector<DataType> &getFieldOnGrid(size_t i) const override {
    return fieldsOnGrids[i];
  }
};

#endif //IC_MULTILEVELFIELD_HPP
