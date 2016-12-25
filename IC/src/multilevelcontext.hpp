#ifndef _MULTILEVELFIELDMANAGER_HPP
#define _MULTILEVELFIELDMANAGER_HPP

#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <cassert>
#include <vector>
#include <memory>

#include "grid.hpp"
#include "signaling.hpp"
#include "field.hpp"
#include "multilevelfield.hpp"

template<typename T>
class ConstraintField;

template<typename T>
class CAMB;

//#include <Eigen/Dense>



template<typename T>
class MultiLevelContextInformation : public Signaling {
private:
  std::vector<std::shared_ptr<Grid<T>>> pGrid;
  std::vector<std::vector<T>> C0s;
  std::vector<T> weights;

protected:
  std::vector<size_t> Ns;
  std::vector<size_t> cumu_Ns;
  size_t Ntot;
  size_t nLevels;
  T simSize;

  void mapIdToLevelId(size_t i, size_t &level, size_t &level_id) {
    level = 0;
    while (level < nLevels && i >= Ns[level]) {
      i -= Ns[level];
      level++;
    }
    if (i >= Ns[level])
      throw std::runtime_error("ID out of range when mapping into underlying fields in mapIdToLevelId");
    level_id = i;
  }

  MultiLevelContextInformation(size_t N) {
    nLevels = 1;
    Ns.push_back(N);
    Ntot = N;
  }


public:

  MultiLevelContextInformation() {
    nLevels = 0;
    Ntot = 0;
  }

  virtual ~MultiLevelContextInformation() {}


  void addLevel(const CAMB<T> &spectrum, T size, size_t nside, const Coordinate<T> & offset={0,0,0}) {
    if(!spectrum.isUsable())
      throw std::runtime_error("Cannot add a grid level until the power spectrum has been specified");

    if(nLevels==0)
      simSize = size;

    auto grid = std::make_shared<Grid<T>>(simSize, nside, size/nside, offset.x, offset.y, offset.z);
    std::vector<T> C0 = spectrum.getPowerSpectrumForGrid(*grid);
    addLevel(C0, grid);

  }

  void addLevel(std::vector<T> C0, std::shared_ptr<Grid<T>> pG) {
    C0s.push_back(std::move(C0));
    if (pGrid.size() == 0) {
      weights.push_back(1.0);
    } else {
      weights.push_back(pow(pG->dx / pGrid[0]->dx, 3.0));
    }
    pGrid.push_back(pG);
    Ns.push_back(pG->size3);
    cumu_Ns.push_back(Ntot);
    Ntot += pG->size3;
    nLevels += 1;
    this->changed();
  }

  void clear() {
    C0s.clear();
    pGrid.clear();
    Ns.clear();
    cumu_Ns.clear();
    Ntot = 0;
    nLevels = 0;
    this->changed();
  }


  size_t getNumCells() const {
    return Ntot;
  }

  size_t getNumLevels() const {
    return nLevels;
  }

  Grid<T> &getGridForLevel(size_t level) {
    return *pGrid[level];
  }

  const Grid<T> &getGridForLevel(size_t level) const {
    return *pGrid[level];
  }

  T getWeightForLevel(size_t level) {
    return weights[level];
  }

  const std::vector<T> &getCovariance(size_t level) const {
    return C0s[level];
  }

  vector<complex<T>> createEmptyFieldForLevel(size_t level) const {
    vector<complex<T>> ar(pGrid[level]->size3);
    return ar;
  }


  auto generateMultilevelFromHighResField(Field<complex<T>, T> &&data) {
    assert(&data.getGrid() == pGrid.back().get());

    // Generate the fields on each level. Fill low-res levels with zeros to start with.
    vector<Field<complex<T>, T>> dataOnLevels;
    for (size_t level = 0; level < pGrid.size(); level++) {
      if (level == pGrid.size() - 1) {
        dataOnLevels.emplace_back(std::move(data));
      } else {
        dataOnLevels.emplace_back(*pGrid[level], false);
      }
    }

    // Now interpolate the high-res level down into lower-res levels
    size_t levelmax = dataOnLevels.size() - 1;
    if (levelmax > 0) {
      assert(dataOnLevels.back().isFourier());
      dataOnLevels.back().toReal();
      for (int level = levelmax - 1; level >= 0; --level) {
        dataOnLevels[level].addFieldFromDifferentGrid(dataOnLevels.back());
      }
    }

    return ConstraintField<std::complex<T>>(*this, std::move(dataOnLevels));
  }

  void forEachLevel(std::function<void(Grid<T> &)> newLevelCallback) {
    for (size_t level = 0; level < nLevels; level++) {
      newLevelCallback(*pGrid[level]);
    }
  }

  void forEachLevel(std::function<void(const Grid<T> &)> newLevelCallback) const {
    for (size_t level = 0; level < nLevels; level++) {
      newLevelCallback(*pGrid[level]);
    }
  }


  void forEachCellOfEachLevel(
    std::function<bool(size_t)> levelCallback,
    std::function<void(size_t, size_t, size_t)> cellCallback,
    bool kspace = true) {


    for (size_t level = 0; level < nLevels; level++) {

      if (!levelCallback(level))
        continue;
      size_t level_base = cumu_Ns[level];

#pragma omp parallel for
      for (size_t i = 0; i < Ns[level]; i++) {
        size_t i_all_levels = level_base + i;
        cellCallback(level, i, i_all_levels);
      }
    }

  }

  void forEachCellOfEachLevel(
    std::function<void(size_t, size_t, size_t, std::vector<std::complex<T>> &)> cellCallback,
    bool kspace = true) {

    forEachCellOfEachLevel([](size_t i) { return true; }, cellCallback);
  }


  std::complex<T> accumulateOverEachCellOfEachLevel(
    std::function<bool(size_t)> newLevelCallback,
    std::function<std::complex<T>(size_t, size_t, size_t)> getCellContribution) {

    T res_real(0), res_imag(0);


    for (size_t level = 0; level < nLevels; level++) {

      if (!newLevelCallback(level))
        continue;

      size_t level_base = cumu_Ns[level];

#pragma omp parallel for reduction(+:res_real,res_imag)
      for (size_t i = 0; i < Ns[level]; i++) {
        size_t i_all_levels = level_base + i;
        std::complex<T> res = getCellContribution(level, i, i_all_levels);

        // accumulate separately - OMP doesn't support complex number reduction :-(
        res_real += std::real(res);
        res_imag += std::imag(res);
      }
    }

    return std::complex<T>(res_real, res_imag);
  }


  void copyContextWithIntermediateResolutionGrids(MultiLevelContextInformation<T> &newStack, size_t base_factor = 2,
                                                  size_t extra_lores = 1) {
    /* Copy this MultiLevelContextInformation, but insert intermediate virtual grids such that
     * there is a full stack increasing in the specified power.
     *
     * E.g. if there is a 256^3 and a 1024^3 grid stack, with default parameters this will return a
     * 128^3, 256^3, 512^3 and 1024^3 grid stack.
     *
     * The first two will be based on the 256^3 literal grid, and the second two on the 1024^3 literal grid.
     *
     * @param newStack     the MultiLevelContextInformation into which the new stack will be placed. Any
     *                     existing grids in the stack will be removed.
     * @param base_factor  grids will be downgraded by factors of base_factor^N where N is an integer
     * @param extra_lores  number of additional grids *below* the base level to add
    */
    newStack.clear();

    for (size_t level = 0; level < nLevels; ++level) {
      size_t neff = size_t(round(pGrid[0]->dx / pGrid[level]->dx)) * pGrid[0]->size;
      if (level > 0) {
        size_t factor = base_factor;
        while (pGrid[level]->dx * factor * 1.001 < pGrid[level - 1]->dx) {
          cerr << "Adding virtual grid with effective resolution " << neff / factor << endl;
          auto vGrid = std::make_shared<SubSampleGrid<T>>(pGrid[level], factor);
          newStack.addLevel(vector<T>(), vGrid);
          factor *= base_factor;
        }
      } else {
        size_t factor = base_factor;
        for (size_t i = 0; i < extra_lores; ++i) {
          cerr << "Adding virtual grid with effective resolution " << neff / factor << endl;
          auto vGrid = std::make_shared<SubSampleGrid<T>>(pGrid[level], factor);
          newStack.addLevel(vector<T>(), vGrid);
          factor *= base_factor;
        }
      }

      cerr << "Adding real grid with resolution " << neff << endl;
      newStack.addLevel(C0s[level], pGrid[level]);

    }

  }

};

#endif
