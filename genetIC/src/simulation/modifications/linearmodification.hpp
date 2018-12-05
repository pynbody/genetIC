#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <cstddef>
#include <vector>
#include <tuple>
#include <src/simulation/grid/grid.hpp>
#include <src/simulation/modifications/modification.hpp>
#include <src/simulation/multilevelcontext/multilevelcontext.hpp>
#include <src/simulation/field/field.hpp>

namespace modifications {


  /*! \brief Class defining linear modifications, ie, those that constrain a linear function of the field vector.
    Prototypical example is the density modification, which fixes the average overdensity in some region to be a particular value.

    Each linear modification is defined by a covector, a, and target d, and acts on a field vector, f, such that
    a.f = d, where . indicates the inner product of the covector and field.

    */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class LinearModification : public Modification<DataType, T> {

  public:
    /*! \brief Constructor, specifying multi level context and cosmological parameters

        \param underlying_ - underlying multi-level context object.
        \param cosmology_ - struct containing cosmological data.
    */
    LinearModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                       const cosmology::CosmologicalParameters<T> &cosmology_) :
        Modification<DataType, T>(underlying_, cosmology_) {
      this->order = 1;

      if(this->flaggedCells[this->underlying.getNumLevels() - 1].empty()){
           throw std::runtime_error("Linear modifications use exclusively information from the highest resolution grid "
                                            "but no cells have been flagged on this grid.\n"
                                            "Change flagged cells selection to use linear modifications.");
      } else {
        this->flaggedCellsFinestGrid = this->flaggedCells[this->underlying.getNumLevels() - 1];
      }

    };

    //! Calculates the current value of a.field for the modification defined by covector a on the field.
    T calculateCurrentValue(const fields::MultiLevelField<DataType> &field) override {
      T val = this->getCovector()->innerProduct(field).real();
      return val;
    }

    //! Returns the covector that defines this modification.
    std::shared_ptr<fields::ConstraintField<DataType>> getCovector() {
      if(this->covector==nullptr) {
        this->covector = this->calculateCovectorOnAllLevels();
        this->covector->toFourier();
      }
      return this->covector;
    }

  protected:
    std::shared_ptr<fields::ConstraintField<DataType>> covector; //!< Linear modification can be described as covectors.
    std::vector<size_t> flaggedCellsFinestGrid; //!< Linear modifications only use high-res information and extrapolate from there.

    //! Calculate covector on finest level and generate from it the multi-grid field
    std::shared_ptr<fields::ConstraintField<DataType>> calculateCovectorOnAllLevels() {

      size_t level = this->underlying.getNumLevels() - 1;

      using tools::numerics::operator/=;
      auto highResModif = calculateLocalisationCovector(this->underlying.getGridForLevel(level));

      if (level != 0) {
        highResModif.getDataVector() /= this->underlying.getWeightForLevel(level);
      }

      // Note - implicitly applies to dark matter only (baryon modifications not implemented):
      auto multiLevelConstraint = this->underlying.generateMultilevelFromHighResField(std::move(highResModif));
      for(level=0; level<this->underlying.getNumLevels(); ++level) {
        turnLocalisationCovectorIntoModificationCovector(multiLevelConstraint->getFieldForLevel(level));
      }
      return multiLevelConstraint;
    }

    //! Returns a covector for the specified grid defined such that a.f returns the average of field f over the flagged points on the grid.
    fields::Field<DataType, T> calculateLocalisationCovector(grids::Grid<T> &grid)  {

      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();


      T w = 1.0 / this->flaggedCellsFinestGrid.size();

      for (size_t i = 0; i < grid.size3; ++i) {
        outputData[i] = 0;
      }

      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
        outputData[this->flaggedCellsFinestGrid[i]] += w;
      }

      outputField.toFourier();
      return outputField;
    }

    //! To be overriden to perform any required non-local manipulation of the covector, e.g. taking gradients etc
    virtual void turnLocalisationCovectorIntoModificationCovector(fields::Field<DataType, T> &fieldOnLevel) const = 0;


    //! Obtain centre of region of interest
    /*!
     * Mostly useful for angular momentum modifications. Has not been tested for two years.
     */
    Coordinate<T> getCentre(grids::Grid<T> &grid) {
      return grid.getCentreWrapped(this->flaggedCellsFinestGrid);
    }

  };


  /*!
  \class OverdensityModification
  \brief Constrains the average of the field at the flagged points to be equal to the target.
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class OverdensityModification : public LinearModification<DataType, T> {
  public:

    /*! \brief Constructor from underlying multi-level context and cosmological data.

        \param underlying_ - underlying multi-level context object.
        \param cosmology_ - struct containing cosmological parameters.
    */
    OverdensityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_) :
                                LinearModification<DataType, T>(underlying_, cosmology_) {

    };


    //! In the case of the overdensity modification, the average of the field is what we want to constrain, so this function doesn't need to do anything.
    void turnLocalisationCovectorIntoModificationCovector(fields::Field<DataType, T> &fieldOnLevel) const override {}
  };


  /*!
  \class PotentialModification
  \brief Constrains the average of the Newtonian potential over the flagged points to be equal to the target
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class PotentialModification : public LinearModification<DataType, T> {
  public:
    /*! \brief Constructor from underlying multi-level context and cosmological data.

        \param underlying_ - underlying multi-level context object.
        \param cosmology_ - struct containing cosmological parameters.
    */
    PotentialModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_) : LinearModification<DataType, T>(
        underlying_, cosmology_) {

    };

    //! Converts the density into a potential, but otherwise will just average over the flagged cells
    void turnLocalisationCovectorIntoModificationCovector(fields::Field<DataType, T> &fieldOnLevel) const override {
      cosmology::densityToPotential(fieldOnLevel, this->cosmology);
    }

  };

  /*! \class VelocityModification
      \brief  Constrains the average velocity along a specified direction.
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class VelocityModification : public OverdensityModification<DataType, T> {
  protected:
    int direction; //!< Component of velocity to modify, (0,1,2) <-> (x,y,x).

  public:
    /*! \brief Construct a velocity modification from a specified multi-level context, cosmology, and velocity component

        \param underlying_ - underlying multi-level context object.
        \param cosmology_ - struct containing cosmological data
        \param direction - component of t velocity to modify, (0,1,2) <-> (x,y,x).
    */
    VelocityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_, int direction_) :
                          OverdensityModification<DataType, T>(underlying_, cosmology_), direction(direction_) {};



      //! Converts the specified field from a multi-level constraint covector into a velocity modification
      void turnLocalisationCovectorIntoModificationCovector(fields::Field<DataType, T> &fieldOnLevel) const override {
        using compT = std::complex<T>;
        fieldOnLevel.toFourier(); // probably already in Fourier space, but best to be sure
        const grids::Grid<T> &grid(fieldOnLevel.getGrid());
        complex<T> I(0, 1);
        T scale = zeldovichVelocityToOffsetRatio(this->cosmology);
        const T nyquist = tools::numerics::fourier::getNyquistModeThatMustBeReal(grid) * grid.getFourierKmin();

        auto calcCell =
                [I, scale, nyquist](complex<T> overdensityFieldValue, T kx, T ky, T kz, T k_chosen_direction) -> complex<T> {

          // This lambda evaluates the derivative for the given Fourier-space cell. kx,ky,kz are the
          // input k-space coordinates, and k_chosen_direction must be set to whichever of these specifies
          // the derivative direction

          T k2 = kx*kx+ky*ky+kz*kz;

          if(k_chosen_direction==nyquist || k2==0)
            return complex<T>(0);
          else
            return -scale * overdensityFieldValue * I * k_chosen_direction/k2;
        };


        if(direction == 0)
          fieldOnLevel.forEachFourierCell([calcCell](compT v, T kx, T ky, T kz) { return calcCell(v, kx, ky, kz, kx);});
        else if(direction == 1)
          fieldOnLevel.forEachFourierCell([calcCell](compT v, T kx, T ky, T kz) { return calcCell(v, kx, ky, kz, ky);});
        else if(direction == 2)
          fieldOnLevel.forEachFourierCell([calcCell](compT v, T kx, T ky, T kz) { return calcCell(v, kx, ky, kz, kz);});
        else
          throw std::runtime_error("Unknown velocity direction");
      }
  };

}

#endif

