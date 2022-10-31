#ifndef IC_MULTILEVELFIELD_HPP
#define IC_MULTILEVELFIELD_HPP

#include "src/simulation/multilevelgrid/multilevelgrid.hpp"
#include "src/simulation/filters/filterfamily.hpp"
#include "src/simulation/field/field.hpp"


namespace fields {

  /*!   \class MultiLevelField
        \brief Class to store a field defined across multiple grids.

        The multi-level field relies on a multi-level context to describe the grids on which the fields are 
        defined. This context is not allowed to change after the field object has been initialised. 
  */
  template<typename DataType>
  class MultiLevelField : public std::enable_shared_from_this<MultiLevelField<DataType>> {

  public:
    using T = tools::datatypes::strip_complex<DataType>;
    using ComplexType = tools::datatypes::ensure_complex<DataType>;

  protected:
    const multilevelgrid::MultiLevelGrid<DataType> *multiLevelContext; //!< Pointer to the underlying multi-level context
    bool isCovector; //!< True if the multi-level field is a covector, used to define a modification.

    std::vector<std::shared_ptr<Field<DataType, T>>> fieldsOnLevels; //!< Vector that stores all the fields on the different levels

    //! \brief Which transfer function the field currently has applied
    particle::species transferType;
 
    //!< Check the multi-level context has not been updated since the field was created
    void assertContextConsistent() const {
      assert(multiLevelContext->getNumLevels() == fieldsOnLevels.size());
    }

    /*! \brief Constructor with fields unspecified - only multi-level context which defines the grids.
     * 
     * This means there are no fields defined on the grid levels, which is an inconsistent state
     * for the object; therefore this should only be called by constructors of child classes.
    */ 
    MultiLevelField(const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext,
                    particle::species transfer_type = particle::species::dm) :
      multiLevelContext(&multiLevelContext), transferType(transfer_type) {
      isCovector = false;
    }
    
  public:

    //! Constructor from fields for each level of a specified multi-level context
    MultiLevelField(const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext,
                    const std::vector<std::shared_ptr<Field<DataType, T>>> &fieldsOnLevels,
                    particle::species transfer_type = particle::species::dm) :
      multiLevelContext(&multiLevelContext), fieldsOnLevels(fieldsOnLevels) {
      assertContextConsistent();
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

    virtual ~MultiLevelField() {}
    
    //! Return the transfer function type currently applied to this field
    particle::species getTransferType() const {
      return transferType;
    }

    //! Returns a reference to the multi-level context associated to this multi-level field.
    virtual multilevelgrid::MultiLevelGrid<DataType> &getContext() const {
      return const_cast<multilevelgrid::MultiLevelGrid<DataType> &>(*multiLevelContext);
    }

    //! Returns a constant reference to the field on level i of the multi-level context.
    virtual const Field<DataType, T> &getFieldForLevel(size_t i) const {
      assertContextConsistent();
      assert(i < fieldsOnLevels.size());
      return *(fieldsOnLevels[i]);
    }

    //! Returns a constant reference to the field on the specified grid.
    virtual const Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) const {
      assertContextConsistent();
      return (const_cast<MultiLevelField<DataType> *>(this))->getFieldForGrid(grid);
    };

    //! Returns a reference to the field on the specified grid.
    virtual Field<DataType, T> &getFieldForGrid(const grids::Grid<T> &grid) {
      assertContextConsistent();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        if (grid.isProxyFor(&multiLevelContext->getGridForLevel(i)))
          return getFieldForLevel(i);
      }
      throw (std::runtime_error("Cannot find a field for the specified grid"));
    };

    //! Returns a reference to the field on level i of the multi-level context (can edit field)
    virtual Field<DataType, T> &getFieldForLevel(size_t i) {
      assertContextConsistent();
      return *(fieldsOnLevels[i]);
    }

    //! Returns the number of levels in this field
    size_t getNumLevels() const {
      return multiLevelContext->getNumLevels();
    }

    //! Checks whether a field is defined on the specified level.
    bool hasFieldForLevel(size_t i) const {
      assertContextConsistent();
      return this->getFieldForLevel(i).getDataVector().size() > 0;
    }

    //! Converts the fields on each level to real space, if they are not already.
    void toReal() {
      assertContextConsistent();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        getFieldForLevel(i).toReal();
    }

    //! Converts the fields on each level to Fourier space, if they are not already.
    void toFourier() {
      assertContextConsistent();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        getFieldForLevel(i).toFourier();
    }

    //! Returns true if the specified field has the same multi-level context as this one.
    bool isCompatible(const MultiLevelField<DataType> &other) const {
      return other.multiLevelContext == multiLevelContext;
    }

    //! Returns true if the field is defined in real space on all levels.
    bool isRealOnAllLevels() const {
      assertContextConsistent();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        if (getFieldForLevel(i).isFourier()) return false;
      return true;
    }

    //! Returns true if the field is defined in Fourier space on all levels.
    bool isFourierOnAllLevels() const {
      assertContextConsistent();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i)
        if (!getFieldForLevel(i).isFourier()) return false;
      return true;
    }

    //! Returns suitable filters for combining information on different grid levels
    filters::FilterFamily<T> getFilters() const {
      assertContextConsistent();
      return filters::FilterFamily<T>(*multiLevelContext);
    }

    //! Adds the specified multi-level field to this one.
    void operator+=(const MultiLevelField<DataType> &other) {
      addScaled(other, 1.0);
    }

    //! Divides the field by the specified ratio.
    void operator/=(DataType ratio) {
      using namespace tools::numerics;

      for (size_t level = 0; level < getNumLevels(); level++) {
        auto &data = getFieldForLevel(level).getDataVector();
        data /= ratio;
      }
    }

    //! Multiplied the field by the specified ratio.
    void operator*=(DataType ratio) {
      using namespace tools::numerics;

      for (size_t level = 0; level < getNumLevels(); level++) {
        auto &data = getFieldForLevel(level).getDataVector();
        data *= ratio;
      }
    }

    //! Flips the sign of the field.
    void reverse() {
      for (size_t level = 0; level < getNumLevels(); ++level) {
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
      assertContextConsistent();
      assert (isCompatible(other));
      if (other.isFourierOnAllLevels()) {
        toFourier();

        for (size_t level = 0; level < getNumLevels(); level++) {
          if (hasFieldForLevel(level) && other.hasFieldForLevel(level)) {
            Field<DataType> &fieldThis = getFieldForLevel(level);
            const Field<DataType> &fieldOther = other.getFieldForLevel(level);
            T kMin = fieldThis.getGrid().getFourierKmin();
            fieldThis.forEachFourierCellInt([&fieldOther, kMin, scale]
                                              (ComplexType currentVal, int kx, int ky, int kz) {
              return currentVal + scale * fieldOther.getFourierCoefficient(kx, ky, kz);
            });
          }
        }
      } else if (other.isRealOnAllLevels()) {
        toReal();
        for (size_t level = 0; level < getNumLevels(); level++) {
          Field<DataType> &fieldThis = getFieldForLevel(level);
          const Field<DataType> &fieldOther = other.getFieldForLevel(level);

          fieldThis.addScaled(fieldOther, scale);
        }
      } else {
        throw std::runtime_error("Expected the other field to be in real/fourier space, not a mix of these.");
      }
    }

    //! \brief Copy across data from another field.
    //! Note that this is different to operator=, which copies all attributes of the object.
    //! Here, we copy only the field data, without acquiring any of its other properties (such as
    //! transferType, which would change the transfer function the field uses).
    void copyData(const MultiLevelField<DataType> &other) {
      assertContextConsistent();
      assert (isCompatible(other));
      assert(other.isFourierOnAllLevels());
      toFourier();
      for (size_t level = 0; level < getNumLevels(); level++) {
        if (hasFieldForLevel(level) && other.hasFieldForLevel(level)) {
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
      applyTransferRatio(other.getTransferType());
    }
    
    //! \brief Takes the inner product between two fields.
    ComplexType innerProduct(const MultiLevelField<DataType> &other) const {

      /* To explain what happens below:
       * The inner product is defined as the operation between a covector a and a vector b : a * b elementwise.
       * Because of the way that covectors and vectors are defined in our approach, this remains true even when the
       * elements are defined on different grids, i.e. the relative cell size ratios and cross-terms between
       * low and high resolution parts are all absorbed into the definition of a covector.
       *
       */
      assertContextConsistent();
      assert(isCompatible(other));
      assert(other.getTransferType()==this->transferType);
      if (!isCovector)
        throw (std::runtime_error(
          "The inner product can only be taken if one of the fields is a covector"));

      assert(isFourierOnAllLevels() && other.isFourierOnAllLevels());
      // To take inner product with correct filters, we must have the fields in fourier space

      if(other.isCovector) {
        // TODO: Potentially, optimise so that conversion is done 'on the fly'
        MultiLevelField<DataType> otherAsVector(other);
        otherAsVector.convertToVector();
        otherAsVector.toFourier();
        return this->innerProduct(otherAsVector);
      }

      const Field<DataType> *pFieldThis, *pFieldOther;
      const std::vector<DataType> *pFieldDataThis;

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

    //! Applies the specified filters to this field
    void applyFilters(const filters::FilterFamilyBase<T> & filters) {
      assertContextConsistent();
      for (size_t level = 0; level < getNumLevels(); ++level) {
        if (hasFieldForLevel(level)) {
          getFieldForLevel(level).applyFilter(filters.getFilterForLevel(level));
        }
      }
    }

    /*! \brief Applies the default filters to this field
     *
     * This has the effect of leaving high-k wavemodes on the high-resolution grids, and low-k wavemodes on the
     * low-resolution grids.
     */
    void applyFilters() {
      applyFilters(this->getFilters());
    }

    /*! \brief Converts the field into a covector, using the covariance matrix associated to the field.

    */
    void convertToCovector() {
      assertContextConsistent();
      assert(!isCovector);
      applyMetric(true);
      isCovector = true;
    }

    //! Converts the field back to a vector if it has been converted to or defined as a covector field.
    void convertToVector() {
      assertContextConsistent();
      assert(isCovector);
      applyMetric(false);
      isCovector = false;
    }

  protected:
    /*! \brief Applies the metric to convert a vector to a covector or vice versa
     *
     * At zeroth order, the inverse metric and forward metric are identical. However, when converting *to* a vector,
     * we assume all high frequency information is fully within the zoom window. This cannot be consistently assumed
     * when making the opposite transformation -- if we are trying to calculate a chi^2, for example, the overdensity
     * field on the low resolution grid does contain information above the filter frequency.
     *
     * Thus, when converting to a covector, a more careful filter-within-window approach is applied based on the
     * analysis given in the notes.
     */
    void applyMetric(bool toCovector = false) {
      auto filters = getFilters();
      toFourier();

      decltype(fieldsOnLevels) newFields;

      for(size_t level=0; level < getNumLevels(); ++level) {
        std::shared_ptr<Field<DataType,T>> result = getFieldForLevel(level).copy();
        newFields.push_back(result);

        auto & f = filters.getFilterForLevel(level);
        if(toCovector && level < getNumLevels()-1) {
          auto window = multiLevelContext->getGridForLevel(level+1).getWindow();
          result->applyFilterInWindow(f, window, true);
          result->applyFilterInWindow(f, window, false);
        } else
          result->applyFilter(f*f);

        for(size_t source_level=0; source_level < getNumLevels(); ++source_level) {
          auto & source_field = getFieldForLevel(source_level);
          T pixel_volume_ratio = multiLevelContext->getWeightForLevel(level) /
                                 multiLevelContext->getWeightForLevel(source_level);

          if(source_level == level) {
            continue; // handled above
          } else {
            // high_from_low term and low_from_high term are both captured by the following expression.
            // Note we do not include the effect of the zoom windowing in this part (even for vector->covector
            // changes) because it is a correction to a correction to a correction... playing with the toy_implementation
            // shows that it has no discernible effect on the measured chi^2 (which is the only place it would enter).
            auto & this_level_filter = filters.getFilterForLevel(level);
            auto & source_level_filter = filters.getFilterForLevel(source_level);
            result->addFieldFromDifferentGridWithFilter(source_field,
              this_level_filter*source_level_filter*sqrt(pixel_volume_ratio));
          }
        }
        result->toFourier();
      }
      fieldsOnLevels = newFields;
      assert(this->isFourierOnAllLevels());
    }

  public:


    /*! \brief Forces 'exact' power spectrum by normalising the white noise field to unit variance.

    This means that effectively only the phase of the field at each Fourier mode is randomised,
    not the amplitude. This can be used to perform simulations that can quickly estimate ensemble parameters
    (see Angulo and Pontzen 2016).
    */
    void enforceUnitVariance() {
      assertContextConsistent();
      toFourier();

      // Currently, enforcing the power spectrum must be done while the field still contains white noise
      // Generalising this would not be hard, but also not necessary.
      assert(this->transferType == particle::species::whitenoise);

      for(size_t i = 0; i < getNumLevels(); ++i) {
        enforcePowerSpectrumOneGrid(getFieldForLevel(i));
      }

    }

    /*! Makes this field suitable for outputting particles of the defined type, by applying a transfer function */
    void applyPowerSpectrumFor(particle::species outputSpecies) {
      applyTransferRatio(this->transferType, outputSpecies);
      this->transferType = outputSpecies;
    }

    //! Divide by power spectrum from an "old" species and multiply by the transfer function for the species of this field
    void applyTransferRatio(particle::species oldSpecies, particle::species newSpecies) {
      assertContextConsistent();

      if (oldSpecies == newSpecies)
        return;

      if (oldSpecies == particle::species::whitenoise ) {
        applyTransfer(newSpecies);
        return;
      }

      toFourier();
      for (size_t i = 0; i < getNumLevels(); ++i) {
        applyTransferRatioOneLevel(oldSpecies, newSpecies, i);
      }

    }

    void applyTransferRatioOneLevel(particle::species oldSpecies, particle::species newSpecies, size_t i) {
      auto &grid = multiLevelContext->getGridForLevel(i);
      if(multiLevelContext->getCovariance(i, oldSpecies)!=nullptr)
        getFieldForLevel(i).applyTransferFunction(*(multiLevelContext->getCovariance(i, oldSpecies)),-0.5);
      if(multiLevelContext->getCovariance(i, newSpecies)!=nullptr)
        getFieldForLevel(i).applyTransferFunction(*(multiLevelContext->getCovariance(i, newSpecies)),0.5);
    }

  protected:

    //! Multiplies the field by the relevant power spectrum in Fourier space.
    void applyTransfer(particle::species ofSpecies) {
      assertContextConsistent();
      toFourier();
      for (size_t i = 0; i < multiLevelContext->getNumLevels(); ++i) {
        auto &grid = multiLevelContext->getGridForLevel(i);

        getFieldForLevel(i).applyTransferFunction(*multiLevelContext->getCovariance(i, ofSpecies), 0.5);

      }
    }


  public:


    //! Returns the value of chi^2, with respect to the relevant covariance matrix.
    T getChi2() const {
      assertContextConsistent();

      if(this->getTransferType()!=particle::species::whitenoise)
        throw std::runtime_error("Cannot calculate the chi^2 after the output field has been generated");
      // Would easily be possible to calculate chi^2 on a single level for an arbitrary transfer function, but for
      // multi-level fields it would presumably require a lot of work (because of the need to implement
      // convertToCovector in this arbitrary case).

      bool returnToReal = false;
      if(!this->isFourierOnAllLevels()) {
        // technically converting to fourier violates const correctness, but we'll put it right below
        const_cast<MultiLevelField<DataType>*>(this)->toFourier();
        returnToReal = true;
      }

      auto self_copy = fields::MultiLevelField<DataType>(*this);
      self_copy.convertToCovector();
      T chi2 = self_copy.innerProduct(*this).real();

      if(returnToReal) {
        // Fix the const violation above
        const_cast<MultiLevelField<DataType>*>(this)->toReal();
      }

      return chi2;
    }


  private:

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
        if(absExistingValue==0.0)
          return std::complex<T>(0.0);
        else
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
    OutputField(const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext,
                particle::species transfer_type)
      : MultiLevelField<DataType>(multiLevelContext, transfer_type) {

      fieldsOnLevelsPopulated = false;
    }

    //! Copy constructor
    OutputField(const OutputField<DataType> &copy) : MultiLevelField<DataType>(copy) {
      fieldsOnLevelsPopulated = copy.fieldsOnLevelsPopulated;
    }



    //! Returns the field on the specified level, populating it if not yet populated.
    Field<DataType, T> &getFieldForLevel(size_t i) override {
      populateFieldsOnLevelsIfRequired();
      return *(this->fieldsOnLevels[i]);
    }

    /* Combine all grids together

    This effectively converts from whitenoise basis into delta basis.
    */

   void combineGrids() {
      auto filters = this->getFilters();
      size_t nlevels = this->getContext().getNumLevels();

      logging::entry() << "Combining information from different levels..." << std::endl;

      for (size_t level = 1; level < nlevels; ++level) {

        // remove the low-frequency information from this level
        this->getFieldForLevel(level).applyFilter(
          filters.getHighPassFilterForLevel(level));

        // replace with the low-frequency information from the level below
        this->getFieldForLevel(level).addFieldFromDifferentGridWithFilter(
          getFieldForLevel(level - 1),
          filters.getLowPassFilterForLevel(level - 1));
      }

      this->getContext().setLevelsAreCombined();
   }

   void combineGridsLowerHigher() {
      auto filters = this->getFilters();
      auto copy = *this;
      size_t nlevels = this->getContext().getNumLevels();

      logging::entry() << "Combining information from different levels..." << std::endl;

      for (size_t level = 1; level < nlevels; ++level) {
        // remove the low-frequency information from this level
        this->getFieldForLevel(level).applyFilter(
          filters.getFilterForLevel(level));

        // replace with the low-frequency information from the level below
        this->getFieldForLevel(level).addFieldFromDifferentGridWithFilter(
          copy.getFieldForLevel(level-1),
          filters.getFilterForLevel(level - 1));
        
        // replace with the high-frequency information from the level above
        if (level < nlevels - 1) {
          this->getFieldForLevel(level).addFieldFromDifferentGridWithFilter(
            copy.getFieldForLevel(level+1),
            filters.getFilterForLevel(level + 1));
        }
      }
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
    \param transferType - transfer function that has been applied to fields on which this covector will act
    */
    ConstraintField(const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext,
                    const std::vector<std::shared_ptr<Field<DataType, T>>> &fieldsOnGrids,
                    particle::species transferType,
                    bool isCovector)
      : MultiLevelField<DataType>(multiLevelContext, std::move(fieldsOnGrids)) {
      this->isCovector = isCovector;
      this->transferType = transferType;
    }



  };
}

#endif
