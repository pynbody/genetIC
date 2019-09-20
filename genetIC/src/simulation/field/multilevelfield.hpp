#ifndef IC_MULTILEVELFIELD_HPP
#define IC_MULTILEVELFIELD_HPP

#include "src/simulation/multilevelcontext/multilevelcontext.hpp"
#include "src/simulation/filters/filterfamily.hpp"
#include "src/simulation/field/field.hpp"


namespace fields {

  /*!   \class MultiLevelField
        \brief Class to manage a field defined across multiple grids.

        The multi-level field it to fields what the multi-level context is to grids. It defines a collection
        of fields on different grids, describing different levels of the field at different resolutions.
  */
  template<typename DataType>
  class MultiLevelField : public std::enable_shared_from_this<MultiLevelField<DataType>> {

  public:
    using T = tools::datatypes::strip_complex<DataType>;
    using ComplexType = tools::datatypes::ensure_complex<DataType>;
  protected:

    const multilevelcontext::MultiLevelContextInformation<DataType> *multiLevelContext; //!< Pointer to the underlying multi-level context
    bool isCovector; //!< True if the multi-level field is a covector, used to define a modification.

    std::vector<std::shared_ptr<Field<DataType, T>>> fieldsOnLevels; //!< Vector that stores all the fields on the different levels



  public:

    //! \brief Variable that stores which transfer function the field should request from the multi-level context.

    particle::species transferType;



    //! Constructor with fields unspecified - only multi-level context.
    MultiLevelField(const multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    particle::species transfer_type = particle::species::dm) :
      multiLevelContext(&multiLevelContext), transferType(transfer_type) {
      isCovector = false;
    }

    //! Constructor with fields and multi-level context. specified.
    MultiLevelField(const multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    const std::vector<std::shared_ptr<Field<DataType, T>>> &fieldsOnGrids,
                    particle::species transfer_type = particle::species::dm) :
      multiLevelContext(&multiLevelContext), fieldsOnLevels(fieldsOnGrids) {
      transferType = transfer_type;
      isCovector = false;
    }

    //! Copy constructor
    MultiLevelField(const MultiLevelField<DataType> &copy) :
      std::enable_shared_from_this<MultiLevelField<DataType>>(), multiLevelContext(&(copy.getContext())) {

      for (size_t level = 0; level < multiLevelContext->getNumLevels(); level++) {
        fieldsOnLevels.push_back(std::make_shared<Field<DataType, T>>(copy.getFieldForLevel(level)));
      }
      transferType = copy.transferType;

      isCovector = copy.isCovector;
    }

    //! Destructor
    virtual ~MultiLevelField() {}

    //! Returns a reference to the multi-level context associated to this multi-level field.
    virtual multilevelcontext::MultiLevelContextInformation<DataType> &getContext() const {
      return const_cast<multilevelcontext::MultiLevelContextInformation<DataType> &>(*multiLevelContext);
    }


    //! Returns a constant reference to the field on level i of the multi-level context (cannot edit field)
    virtual const Field<DataType, T> &getFieldForLevel(size_t i) const {
      assert(i < fieldsOnLevels.size());
      return *(fieldsOnLevels[i]);
    }

    //! Returns a constant reference to the field on the specified grid.
    virtual const Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) const {
      return (const_cast<MultiLevelField<DataType> *>(this))->getFieldForGrid(grid);
    };

    //! Returns a reference to the field on the specified grid.
    virtual Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) {
      // TODO: problematically slow implementation
      // MR: Is this still up to date ? (Oct 2017)
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        if (grid.isProxyFor(&multiLevelContext->getGridForLevel(i)))
          return getFieldForLevel(i);
      }
      throw (std::runtime_error("Cannot find a field for the specified grid"));
    };

    //! Returns a reference to the field on level i of the multi-level context (can edit field)
    virtual Field<DataType, T> &getFieldForLevel(size_t i) {
      return *(fieldsOnLevels[i]);
    }


    //! Returns the number of levels in this field's multi-level-context
    size_t getNumLevels() const {
      return multiLevelContext->getNumLevels();
    }

    //! Checks whether a field is defined on the specified level.
    bool hasFieldOnGrid(size_t i) const {
      return this->getFieldForLevel(i).getDataVector().size() > 0;
    }

    //! Converts the fields on each level to real space, if they are not already.
    void toReal() {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        getFieldForLevel(i).toReal();
    }

    //! Converts the fields on each level to Fourier space, if they are not already.
    void toFourier() {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        getFieldForLevel(i).toFourier();
    }

    //! Returns true if the specified field has the same multi-level context as this one.
    bool isCompatible(const MultiLevelField<DataType> &other) const {
      return other.multiLevelContext == multiLevelContext;
    }

    //! Returns trueif the field is Real space on all levels.
    bool isRealOnAllLevels() const {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        if (getFieldForLevel(i).isFourier()) return false;
      return true;
    }

    //! Returns true if the field is Fourier space on all levels.
    bool isFourierOnAllLevels() const {
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        if (!getFieldForLevel(i).isFourier()) return false;
      return true;
    }

    //! Returns suitable filters for combining information on different grid levels
    filters::FilterFamily<T> getFilters() const {
      return filters::FilterFamily<T>(*multiLevelContext);
    }

    //! Adds the specified multi-level field to this one.
    void operator+=(const MultiLevelField<DataType> &other) {
      assert (isCompatible(other));
      addScaled(other, 1.0);
    }

    //! Divides the field on each level by the specified ratio.
    void operator/=(DataType ratio) {
      using namespace tools::numerics;

      for (size_t level = 0; level < getNumLevels(); level++) {
        auto &data = getFieldForLevel(level).getDataVector();
        data /= ratio;
      }
    }

    //! Flips the sign of the field.
    void reverse() {
      for (size_t level = 0; level < multiLevelContext->getNumLevels(); ++level) {
        auto &field = this->getFieldForLevel(level);
        size_t N = field.getDataVector().size();
        auto &field_data = field.getDataVector();
        for (size_t i = 0; i < N; i++) {
          field_data[i] = -field_data[i];
        }
      }
    }

    //! Add a scaled multilevel field to the current one
    void addScaled(const MultiLevelField &other, DataType scale) {
      assert(other.isFourierOnAllLevels());
      toFourier();

      for (size_t level = 0; level < getNumLevels(); level++) {
        if (hasFieldOnGrid(level) && other.hasFieldOnGrid(level)) {
          Field<DataType> &fieldThis = getFieldForLevel(level);
          const Field<DataType> &fieldOther = other.getFieldForLevel(level);
          T kMin = fieldThis.getGrid().getFourierKmin();
          fieldThis.forEachFourierCellInt([&fieldOther, kMin, scale]
                                            (ComplexType currentVal, int kx, int ky, int kz) {
            T k_value = kMin * sqrt(T(kx * kx) + T(ky * ky) + T(kz * kz));
            return currentVal + scale * fieldOther.getFourierCoefficient(kx, ky, kz);
          });
        }
      }

    }

    //! \brief Copy across data from another field.
    //! Note that this is different to operator=, which copies all attributes of the object.
    //! Here, we copy only the field data, without acquiring any of its other properties (such as
    //! transferType, which would change the transfer function the field uses).
    void copyData(const MultiLevelField<DataType> &other) {
      assert (isCompatible(other));
      assert(other.isFourierOnAllLevels());
      toFourier();
      for (size_t level = 0; level < getNumLevels(); level++) {
        if (hasFieldOnGrid(level) && other.hasFieldOnGrid(level)) {
          Field<DataType> &fieldThis = getFieldForLevel(level);
          const Field<DataType> &fieldOther = other.getFieldForLevel(level);

          if (fieldThis.getGrid() != fieldOther.getGrid()) {
            throw std::runtime_error("Attempting to copy data from incompatible grids");
          }

          auto &dataThis = fieldThis.getDataVector();
          const auto &dataOther = fieldOther.getDataVector();

          assert(dataOther.size() == dataThis.size());
          std::copy(dataOther.begin(), dataOther.end(), dataThis.begin());
        }
      }
    }

    //! \brief Copy across data from another field, making any changes required because of differing particle species
    void copyDataAndUpdateForTransferFunction(const MultiLevelField<DataType> &other) {
      copyData(other);
      applyTransferRatio(other.transferType);
    }

    /*! \brief Takes the inner product between two fields.

     * The inner product is defined as the operation between a covector a and a vector b : a * b elementwise.
     * In our case, a and b are split in Fourier space between high and low k modes using their internal filter.
     * The product reads (a_high * b_high) + (a_low * b_low) + (a_high * b_low) + (a_low * b_high).
     *
     * The key in this approach is that low-k components live on the coarse grid and high-k components live on
     * the fine grid. However, there is a contribution coming from the cross-terms which mixes the different levels.
     * Accounting for this contribution at each level of the multi-level grid is done by multiplying by the filtered
     * b.
     *
     * If the field is in a *recombined* state (i.e. right before output to file, or after output to file), its filters
     * are set up in such a way that the inner product only occurs on the finest zoom level. The innerProduct function
     * respects this. It should therefore return an accurate result e.g. for an overdensity covector which is local
     * to the finest zoom level. But it will be completely wrong for e.g. a velocity covector which is by necessity
     * highly non-local. Fixing this would require implementation of a real-space filter that allows information on
     * the coarse levels to be used outside the zoom window.
     *
     * If the two fields are covectors, an extra weighting is applied to convert one of them to a vector by
     * multiplying by the covariance matrix, i.e. the metric in our space.
     */
    ComplexType innerProduct(const MultiLevelField<DataType> &other) const {

      assert(isCompatible(other));
      if (!isCovector)
        throw (std::runtime_error(
          "The inner product can only be taken if one of the fields is regarded as a covector"));

      assert(isFourierOnAllLevels() && other.isFourierOnAllLevels());
      // To take inner product with correct filters, we must have the fields in fourier space

      bool covariance_weighted = false;

      if(other.isCovector) {
        // TODO: Potentially, optimise so that conversion is done 'on the fly'
        MultiLevelField<DataType> otherAsVector(other);
        otherAsVector.convertToVector();
        return this->innerProduct(otherAsVector);
      }

      const filters::Filter<T> *pFiltOther;
      const grids::Grid<T> *pCurrentGrid;
      const Field<DataType> *pFieldThis, *pFieldOther;
      const std::vector<DataType> *pFieldDataThis;
      const Field<DataType> *pCov;

      ComplexType result(0, 0);

      for (size_t level = 0; level < getNumLevels(); ++level) {
        pFieldThis = &(this->getFieldForLevel(level));
        pFieldDataThis = &(this->getFieldForLevel(level).getDataVector());
        pFieldOther = &(other.getFieldForLevel(level));
        if (pFieldOther != nullptr && pFieldDataThis->size() > 0) {
          result += pFieldThis->accumulateForEachFourierCell(
            [&](tools::datatypes::ensure_complex<DataType> thisFieldVal,
                int kx, int ky, int kz) {
              auto otherFieldVal = pFieldOther->getFourierCoefficient(kx, ky, kz);
              return std::real(std::conj(thisFieldVal) * otherFieldVal);
            });
        }
      }
      return result;
    }

    //! Applies filters on all levels.
    void applyFilters(const filters::FilterFamilyBase<T> & filters) {
      for (size_t level = 0; level < getNumLevels(); ++level) {
        if (hasFieldOnGrid(level)) {
          getFieldForLevel(level).applyFilter(filters.getFilterForLevel(level));
        }
      }
    }

    void applyFilters() {
      applyFilters(this->getFilters());
    }

    /*! \brief Converts the field into a covector, using the covariance matrix associated to the field.

        For this to work, the transferType needs to have been specified for the field (defaulting to 0, dark matter).
        Can be either dark matter (0) or baryonic (1).
    */
    void convertToCovector() {
      assert(!isCovector);

      assert(false); // not implemented!

      /*
      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);

        divideByCovarianceOneGrid(getFieldForLevel(i),
                                  *multiLevelContext->getCovariance(i, this->transferType),
                                  grid,
                                  multiLevelContext->getWeightForLevel(i));

      }
       */
      isCovector = true;
    }

    //! Converts the field back to a vector if it has been converted to or defined as a covector field.
    void convertToVector() {
      using tools::numerics::operator*=;
      assert(isCovector);
      toFourier();


      auto filters = this->getFilters();

      decltype(this->fieldsOnLevels) newFields;

      // Add cross-talk terms
      for(size_t level=0; level<getNumLevels(); ++level) {
        std::shared_ptr<fields::Field<DataType,T>> result = getFieldForLevel(level).copy();
        newFields.push_back(result);

        auto & f = filters.getFilterForLevel(level);
        result->applyFilter(f*f);

        for(size_t source_level=0; source_level < getNumLevels(); ++source_level) {
          auto & source_field = this->getFieldForLevel(source_level);
          T pixel_volume_ratio = multiLevelContext->getWeightForLevel(level)/
                                 multiLevelContext->getWeightForLevel(source_level);

          if(source_level == level) {
            continue; // handled above
          } else if(source_level < level) {
            // high_from_low term
            auto & hpf = filters.getFilterForLevel(level);
            auto & lpf = filters.getFilterForLevel(source_level);
            result->addFieldFromDifferentGridWithFilter(source_field, hpf*lpf*sqrt(pixel_volume_ratio));

          } else if(source_level > level) {
            // low_from_high term
            auto & hpf = filters.getFilterForLevel(source_level);
            auto & lpf = filters.getFilterForLevel(level);
            result->addFieldFromDifferentGridWithFilter(source_field, hpf*lpf*sqrt(pixel_volume_ratio));
          }
        }

      }

      this->fieldsOnLevels = newFields;


      // TODO add cross-talk terms
      isCovector = false;
    }



    /*! \brief Forces 'exact' power spectrum by normalising the white noise field to unit variance.

    This means that effectively only the phase of the field at each Fourier mode is randomised,
    not the amplitude. This can be used to perform simulations that can quickly estimate ensemble parameters
    (see Angulo and Pontzen 2016).
    */
    void enforceExactPowerSpectrum() {
      toFourier();

      // Currently, enforcing the power spectrum must be done while the field still contains white noise
      // Generalising this would not be hard, but also not necessary.
      assert(this->transferType == particle::species::whitenoise);

      for(size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        enforcePowerSpectrumOneGrid(getFieldForLevel(i));
      }

    }

    /*! Makes this field suitable for outputting particles of the defined type, by applying a transfer function */
    void applyPowerSpectrumFor(particle::species outputSpecies) {
      applyTransferRatio(this->transferType, outputSpecies);
    }


  protected:
    //! Divide by power spectrum from an "old" species and multiply by the transfer function for the species of this field
    void applyTransferRatio(particle::species oldSpecies, particle::species newSpecies) {

      if (oldSpecies == newSpecies)
        return;

      if (oldSpecies == particle::species::whitenoise ) {
        applyTransfer(newSpecies);
        this->transferType = newSpecies;
        return;
      }

      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);
        applyTransferRatioOneGrid(getFieldForLevel(i),
                                  *multiLevelContext->getCovariance(i, oldSpecies),
                                  *multiLevelContext->getCovariance(i, newSpecies),
                                  grid);
      }
      this->transferType = newSpecies;
    }

    //! Multiplies the field by the relevant power spectrum in Fourier space.
    void applyTransfer(particle::species ofSpecies) {
      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);

        applySpectrumOneGrid(getFieldForLevel(i),
                             *multiLevelContext->getCovariance(i, ofSpecies),
                             grid);

      }
    }

    //! Divide by one power spectrum and multiply by another for a single level.
    void applyTransferRatioOneGrid(Field<DataType> &field,
                                   const Field<DataType> &spectrum1,
                                   const Field<DataType> &spectrum2,
                                   const grids::Grid<T> &grid) {
      field.forEachFourierCellInt([&grid, &field, &spectrum1, &spectrum2]
                                    (complex<T> existingValue, int kx, int ky, int kz) {
        T sqrt_spec1 = sqrt(spectrum1.getFourierCoefficient(kx, ky, kz).real());
        T sqrt_spec2 = sqrt(spectrum2.getFourierCoefficient(kx, ky, kz).real());
        complex<T> new_val;
        if (sqrt_spec1 == 0.0) {
          assert(sqrt_spec2 == 0.0);
          //This can only happen if the power spectrum is zero, ie k = 0, in which case,
          //it will be zero for baryons too, so just return 0:
          new_val = 0.0;
        } else {
          new_val = existingValue * sqrt_spec2 / sqrt_spec1;
        }
        return new_val;
      });
    }

  public:


    //! Returns the value of chi^2, with respect to the relevant covariance matrix.
    T getChi2() {

      this->toFourier();
      auto self_copy = fields::MultiLevelField<DataType>(*this);

      // TODO: Reinstate
      // self_copy.convertToCovector();
      T chi2 = 0; // self_copy.innerProduct(*this).real();
      return chi2;
    }


  private:
    //! \brief Applies the power spectrum to a single level of the multi level field.
    /*!
    \param field - field to apply to
    \param spectrum - covariance matrix to use
    \param grid - grid on which field is defined.
    */
    void applySpectrumOneGrid(Field<DataType> &field,
                              const Field<DataType> &spectrum,
                              const grids::Grid<T> &grid) {

      field.forEachFourierCellInt([&grid, &field, &spectrum]
                                    (complex<T> existingValue, int kx, int ky, int kz) {
        T sqrt_spec = sqrt(spectrum.getFourierCoefficient(kx, ky, kz).real());
        return existingValue * sqrt_spec;
      });
    }

    //! \brief Invert applySpectrumOneGrid
    /*!
    \param field - field to apply to
    \param spectrum - covariance matrix to use
    \param grid - grid on which field is defined.
    */
    void applyInverseSpectrumOneGrid(Field<DataType> &field,
                                     const Field<DataType> &spectrum,
                                     const grids::Grid<T> &grid) {

      field.forEachFourierCellInt([&grid, &field, &spectrum]
                                    (complex<T> existingValue, int kx, int ky, int kz) {
        T sqrt_spec = sqrt(spectrum.getFourierCoefficient(kx, ky, kz).real());
        if (sqrt_spec == 0.0) {
          return existingValue * 0.0;
        } else {
          return existingValue / sqrt_spec;
        }
      });
    }

    //! \brief Multiplies one level of the multi level field by the relevant covariance matrix.
    /*!
    \param field - field to apply to
    \param spectrum - covariance matrix to use
    \param grid - grid on which field is defined.
    \param weight - extra factor to multiply by along with the covariance.
    */
    void multiplyByCovarianceOneGrid(Field<DataType> &field,
                                     const Field<DataType> &spectrum,
                                     const grids::Grid<T> &grid,
                                     T weight) {

      field.forEachFourierCellInt([weight, &grid, &field, &spectrum]
                                    (complex<T> existingValue, int kx, int ky, int kz) {
        T spec = spectrum.getFourierCoefficient(kx, ky, kz).real() * weight;
        return existingValue * spec;
      });
    }

    //! \brief Divides one level by the relevant covariance matrix.
    /*!
    \param field - field to apply to
    \param spectrum - covariance matrix to use
    \param grid - grid on which field is defined.
    \param weight - extra factor to multiply by along with the covariance.
    */
    void divideByCovarianceOneGrid(Field<DataType> &field,
                                   const Field<DataType> &spectrum,
                                   const grids::Grid<T> &grid,
                                   T weight) {

      field.forEachFourierCellInt([weight, &grid, &field, &spectrum]
                                    (complex<T> existingValue, int kx, int ky, int kz) {
        T spec = spectrum.getFourierCoefficient(kx, ky, kz).real() * weight;
        if (spec == 0) {
          return complex<DataType>(0, 0);
        } else {
          return existingValue / (spec);
        }
      });
    }

    //! Apply 'exact' power spectrum on a single grid:
    /*!
    \param field - field to apply to
    \param spectrum - covariance matrix to use
    \param grid - grid on which field is defined.
    */
    void enforcePowerSpectrumOneGrid(Field<DataType> &field) {
      const grids::Grid<T> &grid = field.getGrid();
      T white_noise_norm = sqrt(T(grid.size3));

      field.forEachFourierCellInt([white_noise_norm, &grid, &field]
                                    (complex<T> existingValue, int kx, int ky, int kz) {
        T absExistingValue = abs(existingValue);
        return existingValue/absExistingValue;
      });
    }

  };

  /*! \class OutputField
      \brief Main type of field used to store overdensities for generating output data.

      This is a type of multi-level field. The main difference to the base class is that is
      stores an integer describing its output state, and has the ability to create the fields on
      each level (zeroed out) if these are requested before they have been setup.
  */

  template<typename DataType>
  class OutputField : public MultiLevelField<DataType> {

  protected:
    using T = typename MultiLevelField<DataType>::T;
    typedef enum {
      PRE_SEPARATION,
      SEPARATED,
      RECOMBINED
    } t_output_state;

    t_output_state outputState; //!< Current output state of the OutputField
    bool fieldsOnLevelsPopulated; //!< True if already populated the fields on all levels.


    //! Populates the fields on each level with all zeros
    void populateFieldsOnLevels() {
      this->multiLevelContext->forEachLevel([this](const grids::Grid<T> &g) {
        this->fieldsOnLevels.emplace_back(std::make_shared<Field<DataType, T>>(g));
      });
      fieldsOnLevelsPopulated = true;
    }

    //! Checks whether fields on each level are populated, and populate them if not.
    void populateFieldsOnLevelsIfRequired() {
      if (!fieldsOnLevelsPopulated)
        populateFieldsOnLevels();
    }

  public:
    //!\param Constructor with a given multi-level context and transfer type.
    /*!
    \param multiLevelContext - multiLevel context to define the field on
    \param transfer_type - specifies the transfer function to use for this field
    */
    OutputField(const multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                particle::species transfer_type)
      : MultiLevelField<DataType>(multiLevelContext, transfer_type) {
      outputState = PRE_SEPARATION;
      fieldsOnLevelsPopulated = false;
    }

    //! Copy constructor
    OutputField(const OutputField<DataType> &copy) : MultiLevelField<DataType>(copy) {
      outputState = copy.outputState;
      fieldsOnLevelsPopulated = copy.fieldsOnLevelsPopulated;
    }


    //! Sets the internal state to RECOMBINED
    void setStateRecombined() {
      // TODO - this can probably be removed
      outputState = RECOMBINED;
    }

    //! Returns the field on the specified level, populating it if not yet populated.
    Field<DataType, T> &getFieldForLevel(size_t i) override {
      populateFieldsOnLevelsIfRequired();
      return *(this->fieldsOnLevels[i]);
    }


  };


  /*! \class ConstraintField
    \brief Fields used for defining constraints.

    Note, ConstraintFields are naturally covectors because they map a field onto a single number. For example;
    a constraint might map the density vector to the mean density in some region.
  */
  template<typename DataType>
  class ConstraintField : public MultiLevelField<DataType> {

  protected:
    using T = typename MultiLevelField<DataType>::T;


  public:
    //! \brief Constructor
    /*!
    \param multiLevelContext - current multi-level context
    \param fieldsOnGrids - fields that define the constraint field on each level
    \param transfer_type - 0 for dark matter, 1 for baryons
    */
    ConstraintField(const multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext,
                    const std::vector<std::shared_ptr<Field<DataType, T>>> &fieldsOnGrids)
      : MultiLevelField<DataType>(multiLevelContext, std::move(fieldsOnGrids)) {
      this->isCovector = true;
    }



  };
}

#endif
