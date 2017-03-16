//
// Created by Andrew Pontzen on 15/06/15.
//

#ifndef IC_CONSTRAINT_HPP
#define IC_CONSTRAINT_HPP

#include <string>
#include <complex>
#include <vector>
#include <map>
#include "src/simulation/grid.hpp"
#include "src/cosmology/parameters.hpp"

namespace constraints {

  /** This class takes responsibility for calculating constraint covectors across a single grid */
  template<typename DataType, typename FloatType=strip_complex<DataType>>
  class ConstraintCalculator {
  protected:
    const grids::Grid<FloatType> &grid;
    fields::Field<DataType, FloatType> outputField;
    std::vector<FloatType> &outputData;

    std::vector<size_t> particleArray;
    const cosmology::CosmologicalParameters<FloatType> &cosmology;
    FloatType x0 = 0., y0 = 0., z0 = 0.; //centre coordinates needed for angmom; set in getCentre()

    void centralDifference4thOrder(long index, int direc, FloatType x0, FloatType y0, FloatType z0) {

      FloatType xp = 0., yp = 0., zp = 0.;
      std::tie(xp, yp, zp) = grid.getCellCentroid(index);

      xp = grid.getWrappedDelta(xp, x0);
      yp = grid.getWrappedDelta(yp, y0);
      zp = grid.getWrappedDelta(zp, z0);

      FloatType c[3] = {0, 0, 0};
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
        FloatType rp = std::sqrt((xp * xp) + (yp * yp) + (zp * zp));
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
        long ind_p1, ind_m1, ind_p2, ind_m2;
        int step1[3] = {0, 0, 0};
        int neg_step1[3] = {0, 0, 0};
        step1[di] = 1;
        neg_step1[di] = -1;

        // N.B. can't wrap - might be on subgrid -> we do it anyway! this will likely break with zoom-in TODO
        //ind_m1=grid.findNextIndNoWrap(index, neg_step1);
        //ind_p1=grid.findNextIndNoWrap(index, step1);
        //ind_m2=grid.findNextIndNoWrap(ind_m1, neg_step1);
        //ind_p2=grid.findNextIndNoWrap(ind_p1, step1);

        ind_m1 = grid.getIndexFromIndexAndStep(index, neg_step1);
        ind_p1 = grid.getIndexFromIndexAndStep(index, step1);
        ind_m2 = grid.getIndexFromIndexAndStep(ind_m1, neg_step1);
        ind_p2 = grid.getIndexFromIndexAndStep(ind_p1, step1);

        FloatType a = -1. / 12. / grid.dx, b = 2. / 3. / grid.dx;  //the signs here so that L ~ - Nabla Phi

        outputData[ind_m2] += (c[di] * a);
        outputData[ind_m1] += (c[di] * b);
        outputData[ind_p1] += (-c[di] * b);
        outputData[ind_p2] += (-c[di] * a);
      }

    }

  public:
    ConstraintCalculator(const grids::Grid<FloatType> &grid,
                         const cosmology::CosmologicalParameters<FloatType> &cosmology)
      : grid(grid), outputField(grid, false),
        outputData(outputField.getDataVector()),
        cosmology(cosmology) {
      grid.getFlaggedCells(particleArray);
      getCentre();
    }

    void getCentre() {

      FloatType xa, ya, za, xb, yb, zb, x0 = 0., y0 = 0., z0 = 0.;

      std::vector<size_t> particleArray;
      grid.getFlaggedCells(particleArray);

      std::tie(xa, ya, za) = grid.getCellCentroid(particleArray[0]);

      for (size_t i = 0; i < particleArray.size(); i++) {
        std::tie(xb, yb, zb) = grid.getCellCentroid(particleArray[i]);

        x0 += grid.getWrappedDelta(xb, xa);
        y0 += grid.getWrappedDelta(yb, ya);
        z0 += grid.getWrappedDelta(zb, za);
      }

      x0 /= particleArray.size();
      y0 /= particleArray.size();
      z0 /= particleArray.size();
      x0 += xa;
      y0 += ya;
      z0 += za;

      this->x0 = x0;
      this->y0 = y0;
      this->z0 = z0;

    }

    void overdensity() {

      FloatType w = 1.0 / particleArray.size();

      for (size_t i = 0; i < grid.size3; ++i) {
        outputData[i] = 0;
      }

      for (size_t i = 0; i < particleArray.size(); i++) {
        outputData[particleArray[i]] += w;
      }

      outputField.toFourier();
    }

    void phi() {
      overdensity();
      densityToPotential(outputField, cosmology);
    }

    void angmom(int direction) {

      if (direction < 0 || direction > 2)
        throw std::runtime_error("Angular momentum direction must be 0 (x), 1 (y) or 2 (z)");

      FloatType x0, y0, z0;
      x0 = this->x0;
      y0 = this->y0;
      z0 = this->z0;

      for (size_t i = 0; i < particleArray.size(); i++) {
        centralDifference4thOrder(particleArray[i], direction, x0, y0, z0);
      }

      outputField.toFourier();

      // The constraint as derived is on the potential. By considering
      // unitarity of FT, we can FT the constraint to get the constraint
      // on the density.
      densityToPotential(outputField, cosmology);

    }

    fields::Field<DataType, FloatType> getResult() {
      return std::move(outputField);
    };

  };


  template<typename DataType, typename FloatType=strip_complex<DataType>>
  fields::Field<DataType, FloatType> calcConstraint(std::string name, const grids::Grid<FloatType> &grid,
                                                    const cosmology::CosmologicalParameters<FloatType> &cosmology) {
    typedef ConstraintCalculator<DataType> CC;

    CC calc(grid, cosmology);

    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

    if (name == "overdensity") {
      calc.overdensity();
    } else if (name == "phi") {
      calc.phi();
    } else if (name == "lx") {
      calc.angmom(0);
    } else if (name == "ly") {
      calc.angmom(1);
    } else if (name == "lz") {
      calc.angmom(2);
    } else {
      throw std::runtime_error("Unknown constraint type " + name);
    }
    return calc.getResult();

  }
}


#endif //IC_CONSTRAINT_HPP
