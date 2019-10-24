#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <cstddef>
#include <vector>
#include <tuple>
#include <src/simulation/grid/grid.hpp>
#include <src/simulation/modifications/modification.hpp>
#include <src/simulation/multilevelgrid/multilevelgrid.hpp>
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
    LinearModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                       const cosmology::CosmologicalParameters<T> &cosmology_) :
      Modification<DataType, T>(underlying_, cosmology_) {
      this->order = 1;

      if (this->flaggedCells[this->underlying.getNumLevels() - 1].empty()) {
        throw std::runtime_error("Linear modifications use exclusively information from the highest resolution grid "
                                 "but no cells have been flagged on this grid.\n"
                                 "Change flagged cells selection to use linear modifications.");
      } else {
        this->flaggedCellsFinestGrid = this->flaggedCells[this->underlying.getNumLevels() - 1];
      }

    };

    //! Calculates the current value of a.field for the modification defined by covector a on the field.
    T calculateCurrentValue(const fields::MultiLevelField<DataType> &field) override {
      T val = this->getCovector(field.getTransferType())->innerProduct(field).real();
      return val;
    }

    //! Returns the covector that defines this modification.
    std::shared_ptr<fields::ConstraintField<DataType>> getCovector(particle::species forSpecies) {
      auto r = this->calculateCovectorOnAllLevels(forSpecies);
      r->toFourier();
      return r;
    }

  protected:
    std::shared_ptr<fields::ConstraintField<DataType>> covector; //!< Linear modification can be described as covectors.
    std::vector<size_t> flaggedCellsFinestGrid; //!< Linear modifications only use high-res information and extrapolate from there.

    //! Calculate covector on finest level and generate from it the multi-grid field
    std::shared_ptr<fields::ConstraintField<DataType>> calculateCovectorOnAllLevels(particle::species forSpecies) {

      size_t level = this->underlying.getNumLevels() - 1;

      using tools::numerics::operator/=;
      auto highResModif = calculateLocalisationCovector(this->underlying.getGridForLevel(level));

      auto multiLevelCovector = this->underlying.generateMultilevelCovectorFromHiresCovector(std::move(highResModif), forSpecies);

      turnLocalisationCovectorIntoModificationCovector(*multiLevelCovector);
      return multiLevelCovector;
    }

    //! Returns a covector for the specified grid defined such that a.f returns the average of field f over the flagged points on the grid.
    virtual fields::Field<DataType, T> calculateLocalisationCovector(const grids::Grid<T> &grid) {
      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();

      T w = 1.0 / this->flaggedCellsFinestGrid.size();

      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
        outputData[this->flaggedCellsFinestGrid[i]] += w;
      }

      outputField.toFourier();
      return outputField;
    }

    //! To be overriden to perform any required non-local manipulation of the covector, e.g. taking gradients etc
    virtual void
    turnLocalisationCovectorIntoModificationCovector(fields::MultiLevelField<DataType> & field) const = 0;


    //! Obtain centre of region of interest
    /*!
     * Mostly useful for angular momentum modifications. Has not been tested for two years.
     */
    Coordinate<T> getCentre(const grids::Grid<T> &grid) const {
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
    OverdensityModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_) :
      LinearModification<DataType, T>(underlying_, cosmology_) {

    };


    /*! The modifications are normally made to the white noise field, so the covector itself needs to carry the power spectrum
     *
     * If we are calculating this covector for e.g. the DM field rather than the white noise field, this will
     * automatically be a null-op since we will already have set our particle species to particle::species::dm (in
     * calculateCovectorOnAllLevels, above).
     */
    void turnLocalisationCovectorIntoModificationCovector(fields::MultiLevelField<DataType> & field) const override
    {
      field.applyTransferRatio(field.getTransferType(), particle::species::dm);
      // Note that this does NOT update the field's own transferType to dm, which is the correct behaviour
      // e.g. it may be a covector that acts on whitenoise vectors, in which case we pick up a factor of the
      // DM transfer function here (so that the output is a DM overdensity). The transferType of the covector
      // refers to which type of vector it will act on.
    }
  };


  /*!
  \class PotentialModification
  \brief Constrains the average of the Newtonian potential over the flagged points to be equal to the target
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class PotentialModification : public OverdensityModification<DataType, T> {
  public:
    /*! \brief Constructor from underlying multi-level context and cosmological data.

        \param underlying_ - underlying multi-level context object.
        \param cosmology_ - struct containing cosmological parameters.
    */
    PotentialModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_) : OverdensityModification<DataType, T>(
      underlying_, cosmology_) {

    };

    //! Converts the density into a potential, but otherwise will just average over the flagged cells
    void turnLocalisationCovectorIntoModificationCovector(fields::MultiLevelField<DataType> & field) const override {
      // TODO: combine these operations to make more efficient
      OverdensityModification<DataType,T>::turnLocalisationCovectorIntoModificationCovector(field);
      for(size_t level=0; level<field.getNumLevels(); ++level) {
        auto &fieldOnLevel = field.getFieldForLevel(level);
        cosmology::densityToPotential(fieldOnLevel, this->cosmology);
      }
    }

  };

  /*! \class VelocityModification
      \brief  Constrains the average velocity along a specified direction.
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class VelocityModification : public OverdensityModification<DataType, T> {
  private:
    bool derivInRealSpace = true; //!< If true, make the derivation in real space. Otherwise, make the derivation in Fourier space.
  protected:
    int direction; //!< Component of velocity to modify, (0,1,2) <-> (x,y,z).

  public:
    /*! \brief Construct a velocity modification from a specified multi-level context, cosmology, and velocity component

        \param underlying_ - underlying multi-level context object.
        \param cosmology_ - struct containing cosmological data
         \param direction - component of t velocity to modify, (0,1,2) <-> (x,y,z).
    */
    VelocityModification(const multilevelgrid::MultiLevelGrid<DataType> &underlying_,
                         const cosmology::CosmologicalParameters<T> &cosmology_, int direction_) :
      OverdensityModification<DataType, T>(underlying_, cosmology_), direction(direction_) {};

  protected:

    fields::Field<DataType, T> calculateLocalisationCovector(const grids::Grid<T> &grid) override {
      if (!this->derivInRealSpace) {
        return OverdensityModification<DataType,T>::calculateLocalisationCovector(grid);
      } else {
        // Uses a finite difference 4th order stencil to create just a derivative covector
        Coordinate<int> directionVector;
        Coordinate<int> negDirectionVector;

        directionVector[direction] = 1;
        negDirectionVector[direction] = -1;

        fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
        std::vector<DataType> &outputData = outputField.getDataVector();

        T w = 1.0 / this->flaggedCellsFinestGrid.size();

        size_t ind_p1, ind_m1, ind_p2, ind_m2;

        // Coeffs for the finite difference.  The signs here so that result is - Nabla Phi
        T a = -w / 12. / grid.cellSize, b = w * 2. / 3. / grid.cellSize;

        for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
          size_t index = this->flaggedCellsFinestGrid[i];
          ind_m1 = grid.getIndexFromIndexAndStep(index, negDirectionVector);
          ind_p1 = grid.getIndexFromIndexAndStep(index, directionVector);
          ind_m2 = grid.getIndexFromIndexAndStep(ind_m1, negDirectionVector);
          ind_p2 = grid.getIndexFromIndexAndStep(ind_p1, directionVector);
          outputData[ind_m2] += a;
          outputData[ind_m1] += b;
          outputData[ind_p1] -= b;
          outputData[ind_p2] -= a;
        }

        outputField.toFourier();
        return outputField;

      }

    }

    void computeVelocity(fields::MultiLevelField<DataType> & field, int direction) const {
      for(size_t level=0; level<field.getNumLevels(); ++level) {
        auto &fieldOnLevel = field.getFieldForLevel(level);
        using compT = std::complex<T>;
        fieldOnLevel.toFourier(); // probably already in Fourier space, but best to be sure
        const grids::Grid<T> &grid(fieldOnLevel.getGrid());
        complex<T> I(0, 1);
        T scale = zeldovichVelocityToOffsetRatio(this->cosmology);
        const T nyquist = tools::numerics::fourier::getNyquistModeThatMustBeReal(grid) * grid.getFourierKmin();

        auto calcCell =
          [I, scale, nyquist](complex<T> overdensityFieldValue, T kx, T ky, T kz, T k_chosen_direction, bool derivInRealSpace) -> complex<T> {

            // This lambda evaluates the derivative for the given Fourier-space cell. kx,ky,kz are the
            // input k-space coordinates, and k_chosen_direction must be set to whichever of these specifies
            // the derivative direction

            T k2 = kx * kx + ky * ky + kz * kz;

            if (k2 == 0)
              return complex<T>(0);
            else {
              if (!derivInRealSpace) {
                if(k_chosen_direction == nyquist) return complex<T>(0);
                return -scale * overdensityFieldValue * I * k_chosen_direction / k2;
              } else {
                return -scale * overdensityFieldValue / k2;
              }
            }
          };

        if (direction == 0)
          fieldOnLevel.forEachFourierCell(
            [calcCell, this](compT v, T kx, T ky, T kz) { return calcCell(v, kx, ky, kz, kx, this->derivInRealSpace); });
        else if (direction == 1)
          fieldOnLevel.forEachFourierCell(
            [calcCell, this](compT v, T kx, T ky, T kz) { return calcCell(v, kx, ky, kz, ky, this->derivInRealSpace); });
        else if (direction == 2)
          fieldOnLevel.forEachFourierCell(
            [calcCell, this](compT v, T kx, T ky, T kz) { return calcCell(v, kx, ky, kz, kz, this->derivInRealSpace); });
        else
          throw std::runtime_error("Unknown direction");
      }
    }

  public:

    //! Converts (in-place) from a overdensity covector to a velocity covector
    void turnLocalisationCovectorIntoModificationCovector(fields::MultiLevelField<DataType> & field) const override {
      if(this->underlying.getLevelsAreCombined())
        throw std::runtime_error("Cannot create velocity covector for recombined fields (i.e. after output generation is complete)");
      // Note on the above exception: the problem is that we need to use the Fourier filters to calculate the
      // multi-level Poisson solution. Once Fourier modes are combined on zoom grids, there is no simple way to get the potential.
      OverdensityModification<DataType,T>::turnLocalisationCovectorIntoModificationCovector(field);
      this->computeVelocity(field, this->direction);
    }
  };

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class AngMomentumModification : public OverdensityModification<DataType, T> {
  protected:
    int direction;

  public:
    AngMomentumModification(const multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_, int direction_) :
        OverdensityModification<DataType, T>(underlying_, cosmology_), direction(direction_) {
      if (direction_ < 0 || direction_ > 2)
        throw std::runtime_error("Angular momentum direction must be 0 (x), 1 (y) or 2 (z)");

      direction = direction_;

    };

  protected:
    fields::Field<DataType, T> calculateLocalisationCovector(const grids::Grid<T> &grid) override {

      auto x0 = this->getCentre(grid);
      T a = -1. / 12. / grid.cellSize, b = 2. / 3. / grid.cellSize;  // signs here so that L ~ - Nabla Phi

      int dirp1 = (direction + 1) % 3,
          dirp2 = (direction + 2) % 3;

      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();

      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); ++i) {
        size_t index = this->flaggedCellsFinestGrid[i];
        Coordinate<T> xp = grid.getCentroidFromIndex(index);
        Coordinate<T> dx = grid.getWrappedOffset(xp, x0);
        Coordinate<T> rCrossCoeff;
        // Coefficient to compute cross product (this is the "r \cross" part of r \cross v)
        rCrossCoeff[direction] = 0;
        rCrossCoeff[dirp2] =  dx[dirp1];  // lx = ry * vz
        rCrossCoeff[dirp1] = -dx[dirp2];  //              - rz * vy

        // Create gradient + cross product covector by combining gradient operator with cross product
        // to yield -r \cross \nabla · operator (using 4-th order finite-difference)
        for (int dir = 0; dir < 3; ++dir) {
          if (dir != direction) { // Not necessary, but saves one iteration
            size_t ind_p1, ind_m1, ind_p2, ind_m2;
            Coordinate<int> directionVector;
            Coordinate<int> negDirectionVector;

            directionVector[dir] = 1;
            negDirectionVector[dir] = -1;

            ind_m1 = grid.getIndexFromIndexAndStep(index, negDirectionVector);
            ind_m2 = grid.getIndexFromIndexAndStep(ind_m1, negDirectionVector);
            ind_p1 = grid.getIndexFromIndexAndStep(index, directionVector);
            ind_p2 = grid.getIndexFromIndexAndStep(ind_p1, directionVector);
            outputData[ind_m2] += rCrossCoeff[dir] * a;
            outputData[ind_m1] += rCrossCoeff[dir] * b;
            outputData[ind_p1] -= rCrossCoeff[dir] * b;
            outputData[ind_p2] -= rCrossCoeff[dir] * a;
          }
        }
      }

      outputField.toFourier();

      return outputField;
    }

  public:
    void turnLocalisationCovectorIntoModificationCovector (fields::MultiLevelField<DataType> & field) const override {
      OverdensityModification<DataType,T>::turnLocalisationCovectorIntoModificationCovector(field);

      field.toFourier();
      for(size_t level=0; level<field.getNumLevels(); ++level) {
        auto &fieldOnLevel = field.getFieldForLevel(level);

        complex<T> I(0, 1);
        T scale = zeldovichVelocityToOffsetRatio(this->cosmology);

        auto calcCell =
          [I, scale](complex<T> overdensityFieldValue, T kx, T ky, T kz) -> complex<T> {

            T k2 = kx * kx + ky * ky + kz * kz;

            if (k2 == 0)
              return complex<T>(0);
            else
              return -scale * overdensityFieldValue / k2;
        };

        // Missing k^2 factor (density->potential) and scale (potential->velocity)
        fieldOnLevel.forEachFourierCell(calcCell);

      }
    }
  };

}

#endif

