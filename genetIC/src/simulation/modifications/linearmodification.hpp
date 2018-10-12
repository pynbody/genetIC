#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/modifications/modification.hpp>

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
      T val = this->covector->innerProduct(field).real();
      return val;
    }

    std::shared_ptr<fields::ConstraintField<DataType>> getCovector() {
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
      auto highResModif = calculateCovectorOnOneLevel(this->underlying.getGridForLevel(level));

      if (level != 0) {
        highResModif.getDataVector() /= this->underlying.getWeightForLevel(level);
      }
      return this->underlying.generateMultilevelFromHighResField(std::move(highResModif),this->transferType);
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


    /*!
     * Mostly useful for angular momentum modifications. Has not been tested for two years.
     */
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

        // N.B. can't wrap - might be on subgrid -> we do it anyway! this will likely break with zoom-in
        //ind_m1=grid.findNextIndNoWrap(index, neg_step1);
        //ind_p1=grid.findNextIndNoWrap(index, step1);
        //ind_m2=grid.findNextIndNoWrap(ind_m1, neg_step1);
        //ind_p2=grid.findNextIndNoWrap(ind_p1, step1);

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
                            const cosmology::CosmologicalParameters<T> &cosmology_,size_t transfer_type) :
        LinearModification<DataType, T>(underlying_, cosmology_,transfer_type) {
      this->covector = this->calculateCovectorOnAllLevels();
      this->covector->toFourier();
    };

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) override {

      //Want to return the average of the field over only the flagged cells.
      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();


      //1/number of flagged particles:
      T w = 1.0 / this->flaggedCellsFinestGrid.size();

      //Include only flagged cells in the sum, so first zero out the covector:
      for (size_t i = 0; i < grid.size3; ++i) {
        outputData[i] = 0;
      }

      //Then weight every flagged cell by 1/(# of flagged cells):
      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
        outputData[this->flaggedCellsFinestGrid[i]] += w;
      }

      outputField.toFourier();
      return outputField;
    }
  };


  /*!
  \class PotentialModification
  \brief
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class PotentialModification : public LinearModification<DataType, T> {
  public:
    PotentialModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
                          const cosmology::CosmologicalParameters<T> &cosmology_,size_t transfer_type) : LinearModification<DataType, T>(
        underlying_, cosmology_,transfer_type) {
      this->covector = this->calculateCovectorOnAllLevels();
      this->covector->toFourier();
    };

    fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<T> &grid) override {
      fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid, false);
      std::vector<DataType> &outputData = outputField.getDataVector();

      //Weight for all flagged cells is 1/(# of flagged cells), used to take average over
      //flagged cells only:
      T w = 1.0 / this->flaggedCellsFinestGrid.size();

      //Zero out everythin so we only include flagged cells:
      for (size_t i = 0; i < grid.size3; ++i) {
        outputData[i] = 0;
      }

      //Apply wweight to flagged cells:
      for (size_t i = 0; i < this->flaggedCellsFinestGrid.size(); i++) {
        outputData[this->flaggedCellsFinestGrid[i]] += w;
      }

      //Convert overdensity into a Newtonian potential:
      densityToPotential(outputField, this->cosmology);
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
                            const cosmology::CosmologicalParameters<T> &cosmology_, int direction_,size_t transfer_type) :
        LinearModification<DataType, T>(underlying_, cosmology_,transfer_type) {

      if (direction_ < 0 || direction_ > 2)
        throw std::runtime_error("Angular momentum direction must be 0 (x), 1 (y) or 2 (z)");

      direction = direction_;
      this->covector = this->calculateCovectorOnAllLevels();
      this->covector->toFourier();
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

