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


  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class LinearModification : public Modification<DataType, T> {

  public:

    LinearModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                       const cosmology::CosmologicalParameters<T> &cosmology_,size_t transfer_type) :
        Modification<DataType, T>(underlying_, cosmology_,transfer_type) {
      this->order = 1;
      this->flaggedCellsFinestGrid = this->flaggedCells[this->underlying.getNumLevels() - 1];
    };

    T calculateCurrentValue(const fields::MultiLevelField<DataType> &field) override {
      T val = this->getCovector()->innerProduct(field).real();
      return val;
    }

    std::shared_ptr<fields::ConstraintField<DataType>> getCovector() {
      if(this->covector==nullptr) {
        this->covector = this->calculateCovectorOnAllLevels();
        this->covector->toFourier();
      }
      return this->covector;
    }

  protected:
    std::shared_ptr<fields::ConstraintField<DataType>> covector;      /*!< Linear modification can be described as covectors */
    std::vector<size_t> flaggedCellsFinestGrid;       /*!< Linear modifications only use high-res information and extrapolate from there */

    //! Calculate covector on finest level and generate from it the multi-grid field
    //! \param nField - type of field to create. nField = 0 for dark matter (default), 1 for baryons
    std::shared_ptr<fields::ConstraintField<DataType>> calculateCovectorOnAllLevels() {

      size_t level = this->underlying.getNumLevels() - 1;

      using tools::numerics::operator/=;
      auto highResModif = calculateLocalisationCovector(this->underlying.getGridForLevel(level));

      if (level != 0) {
        highResModif.getDataVector() /= this->underlying.getWeightForLevel(level);
      }

      //Note - implicitly applies to dark matter only (baryon modifications not implemented):
      auto multiLevelConstraint = this->underlying.generateMultilevelFromHighResField(std::move(highResModif));
      for(level=0; level<this->underlying.getNumLevels(); ++level) {
        turnLocalisationCovectorIntoModificationCovector(multiLevelConstraint->getFieldForLevel(level));
      }
      return multiLevelConstraint;
    }

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


    void
    centralDifference4thOrder(grids::Grid<T> &grid, std::vector<DataType> &outputData, size_t index, int direc, T x0,
                              T y0, T z0) {

      T xp = 0., yp = 0., zp = 0.;
      std::tie(xp, yp, zp) = grid.getCentroidFromIndex(index);

      xp = grid.getWrappedOffset(xp, x0);
      yp = grid.getWrappedOffset(yp, y0);
      zp = grid.getWrappedOffset(zp, z0);

      T c[3] = {0, 0, 0};
      if (direc == 0) {
        c[2] = yp;
        c[1] = -zp;
      } else if (direc == 1) {
        c[0] = zp;
        c[2] = -xp;
      } else if (direc == 2) {
        c[1] = xp;
        c[0] = -yp;
      } else if (direc == 3) {
        T rp = std::sqrt((xp * xp) + (yp * yp) + (zp * zp));
        if (rp != 0) {
          c[0] = xp / rp;
          c[1] = yp / rp;
          c[2] = zp / rp;
        }
      } // radial velocity

      else {
        throw std::runtime_error("Wrong value for parameter 'direc' in function 'cen_deriv4_alpha'");
      }

      for (int di = 0; di < 3; di++) {
        size_t ind_p1, ind_m1, ind_p2, ind_m2;
        int step1[3] = {0, 0, 0};
        int neg_step1[3] = {0, 0, 0};
        step1[di] = 1;
        neg_step1[di] = -1;

        ind_m1 = grid.getIndexFromIndexAndStep(index, neg_step1);
        ind_p1 = grid.getIndexFromIndexAndStep(index, step1);
        ind_m2 = grid.getIndexFromIndexAndStep(ind_m1, neg_step1);
        ind_p2 = grid.getIndexFromIndexAndStep(ind_p1, step1);

        T a = -1. / 12. / grid.cellSize, b = 2. / 3. / grid.cellSize;  //the signs here so that L ~ - Nabla Phi

        outputData[ind_m2] += (c[di] * a);
        outputData[ind_m1] += (c[di] * b);
        outputData[ind_p1] += (-c[di] * b);
        outputData[ind_p2] += (-c[di] * a);
      }

    }
  };


  /*!
  \class OverdensityModification
  \brief Constrains the average of the field at the flagged points to be equal to the target.
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class OverdensityModification : public LinearModification<DataType, T> {
  public:

    OverdensityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_) :
                                LinearModification<DataType, T>(underlying_, cosmology_) {

    };


    void turnLocalisationCovectorIntoModificationCovector(fields::Field<DataType, T> &fieldOnLevel) const override {}
  };


  /*!
  \class PotentialModification
  \brief
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class PotentialModification : public LinearModification<DataType, T> {
  public:
    PotentialModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_) : LinearModification<DataType, T>(
        underlying_, cosmology_) {

    };

    void turnLocalisationCovectorIntoModificationCovector(fields::Field<DataType, T> &fieldOnLevel) const override {
      cosmology::densityToPotential(fieldOnLevel, this->cosmology);
    }

  };

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class VelocityModification : public OverdensityModification<DataType, T> {
  protected:
    int direction;

  public:
    VelocityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_, int direction_) :
                          OverdensityModification<DataType, T>(underlying_, cosmology_), direction(direction_) {};



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


  //! WARNING : Unfinished and not working implementation TODO
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class AngMomentumModification : public LinearModification<DataType, T> {
  public:
    int direction;

  public:
    AngMomentumModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_, int direction_,size_t transfer_type) :
        LinearModification<DataType, T>(underlying_, cosmology_,transfer_type) {

      if (direction_ < 0 || direction_ > 2)
        throw std::runtime_error("Angular momentum direction must be 0 (x), 1 (y) or 2 (z)");

      direction = direction_;

    };

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid)  {


      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();


      T x0, y0, z0;
      x0 = this->getCentre(grid).x;
      y0 = this->getCentre(grid).y;
      z0 = this->getCentre(grid).z;

      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
        this->centralDifference4thOrder(grid, outputData, this->flaggedCellsFinestGrid[i], this->direction, x0, y0, z0);
      }

      outputField.toFourier();

      // The modification as derived is on the potential. By considering
      // unitarity of FT, we can FT it to obtain the modified density
      densityToPotential(outputField, this->cosmology);

      return outputField;
    }
  };

}

#endif

