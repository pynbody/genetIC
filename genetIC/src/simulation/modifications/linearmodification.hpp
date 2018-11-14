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
                       const cosmology::CosmologicalParameters<T> &cosmology_) :
        Modification<DataType, T>(underlying_, cosmology_) {
      this->order = 1;
      this->flaggedCellsFinestGrid = this->flaggedCells[this->underlying.getNumLevels() - 1];
    };

    T calculateCurrentValue(const fields::MultiLevelField<DataType> &field) override {
      T val = this->getCovector()->innerProduct(field).real();
      return val;
    }

    std::shared_ptr<fields::ConstraintField<DataType>> getCovector() {
      if(this->covector==nullptr)
        this->covector = this->calculateCovectorOnAllLevels();
      return this->covector;
    }

  protected:
    std::shared_ptr<fields::ConstraintField<DataType>> covector;      /*!< Linear modification can be described as covectors */
    std::vector<size_t> flaggedCellsFinestGrid;       /*!< Linear modifications only use high-res information and extrapolate from there */

    //! Calculate covector on finest level and generate from it the multi-grid field
    std::shared_ptr<fields::ConstraintField<DataType>> calculateCovectorOnAllLevels() {

      size_t level = this->underlying.getNumLevels() - 1;

      using tools::numerics::operator/=;
      auto highResModif = calculateCovectorOnOneLevel(this->underlying.getGridForLevel(level));

      if (level != 0) {
        highResModif.getDataVector() /= this->underlying.getWeightForLevel(level);
      }
      return this->underlying.generateMultilevelFromHighResField(std::move(highResModif));
    }

    //! To be overridden with specific implementations of linear properties
    virtual fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) = 0;


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


  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class OverdensityModification : public LinearModification<DataType, T> {
  public:

    OverdensityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_) :
        LinearModification<DataType, T>(underlying_, cosmology_) {

    };

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) override {

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
  };


  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class PotentialModification : public LinearModification<DataType, T> {
  public:
    PotentialModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_) : LinearModification<DataType, T>(
        underlying_, cosmology_) {

    };

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) override {
      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();

      T w = 1.0 / this->flaggedCellsFinestGrid.size();

      for (size_t i = 0; i < grid.size3; ++i) {
        outputData[i] = 0;
      }

      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
        outputData[this->flaggedCellsFinestGrid[i]] += w;
      }

      densityToPotential(outputField, this->cosmology);
      return outputField;

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

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) override {
      fields::Field<DataType, T> outputField = OverdensityModification<DataType, T>::calculateCovectorOnOneLevel(grid);
      outputField.toFourier(); // probably already in Fourier space, but best to be sure
      std::complex<T> I(0,1);
      T scale = cosmology::zeldovichVelocityToOffsetRatio(this->cosmology);

      auto kinv = [](T kx, T ky, T kz) -> T {
          T k = std::sqrt(kx*kx+ky*ky+kz*kz);
          if(k==0)
            return 0;
          else
            return 1./k;
      };

      if(direction==0)
        outputField.forEachFourierCell([I, kinv, scale](std::complex<T> current_value, T kx, T ky, T kz) {

          return scale * current_value * I * kx*kinv(kx,ky,kz);
        });
      else if(direction==1)
        outputField.forEachFourierCell([I, kinv, scale](std::complex<T> current_value, T kx, T ky, T kz) {
          return scale * current_value * I * ky*kinv(kx,ky,kz);
        });
      else if(direction==2)
        outputField.forEachFourierCell([I, kinv, scale](std::complex<T> current_value, T kx, T ky, T kz) {
          return scale * current_value * I * kz*kinv(kx,ky,kz);
        });
      else {
        throw std::runtime_error("Unknown direction passed to velocity modification");
      }
      return outputField;

    }
  };


  //! WARNING : Unfinished and not working implementation TODO
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class AngMomentumModification : public LinearModification<DataType, T> {
  public:
    int direction;

  public:
    AngMomentumModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                            const cosmology::CosmologicalParameters<T> &cosmology_, int direction_) :
        LinearModification<DataType, T>(underlying_, cosmology_) {

      if (direction_ < 0 || direction_ > 2)
        throw std::runtime_error("Angular momentum direction must be 0 (x), 1 (y) or 2 (z)");

      direction = direction_;

    };

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) override {


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

