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

  virtual ~MultiLevelContextInformation() { }



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


  vector<vector<complex<T>>> generateMultilevelFromHighResField(vector<complex<T>> &&data) {
    assert(data.size() == pGrid.back()->size3);
    vector<vector<complex<T>>> dataOnLevels;
    for (size_t level = 0; level < pGrid.size(); level++) {
      if (level == pGrid.size() - 1) {
        dataOnLevels.emplace_back(std::move(data));
      } else {
        dataOnLevels.emplace_back(pGrid[level]->size3, 0.0);
      }
    }
    size_t levelmax = dataOnLevels.size() - 1;
    if (levelmax > 0) {
      fft(dataOnLevels.back(), dataOnLevels.back(), -1);
      for (int level = levelmax - 1; level >= 0; --level) {

        pGrid[level]->addFieldFromDifferentGrid(*pGrid[levelmax], dataOnLevels[levelmax],
                                                dataOnLevels[level]);


        fft(dataOnLevels[level], dataOnLevels[level], 1);


      }
      fft(dataOnLevels.back(), dataOnLevels.back(), 1);


    }

    return dataOnLevels;

  }



  void applyCovarianceToWhiteNoiseField() {


  }

  void zeroLevel(int level) {
    auto &field = pGrid[level]->getField();
    std::fill(field.begin(), field.end(), 0);
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
    std::function<void(size_t, size_t, size_t, std::vector<std::complex<T>> &)> cellCallback,
    bool kspace = true) {


    for (size_t level = 0; level < nLevels; level++) {

      if(!levelCallback(level))
        continue;
      size_t level_base = cumu_Ns[level];
      auto &field = kspace ? pGrid[level]->getFieldFourier() : pGrid[level]->getField();

#pragma omp parallel for
      for (size_t i = 0; i < Ns[level]; i++) {
        size_t i_all_levels = level_base + i;
        cellCallback(level, i, i_all_levels, field);
      }
    }

  }

  void forEachCellOfEachLevel(
    std::function<void(size_t, size_t, size_t, std::vector<std::complex<T>> &)> cellCallback,
    bool kspace = true) {

    forEachCellOfEachLevel([](size_t i) { return true;}, cellCallback);
  }


  std::complex<T> accumulateOverEachCellOfEachLevel(
    std::function<bool(size_t)> newLevelCallback,
    std::function<std::complex<T>(size_t, size_t, size_t)> getCellContribution) {

    T res_real(0), res_imag(0);


    for (size_t level = 0; level < nLevels; level++) {

      if(!newLevelCallback(level))
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

  /*
  std::complex<T> accumulateOverEachCellOfEachLevel(
    std::function<std::complex<T>(size_t, size_t, size_t, T)> getCellContribution) {
    T weight;
    auto newLevelCallback = [&weight, this](size_t comp) {
      weight = weights[comp];
    };
    auto getCellContributionWrapper = [&weight, &getCellContribution](size_t component, size_t i, size_t cumu_i) {
      return getCellContribution(component, i, cumu_i, weight);
    };
    return accumulateOverEachCellOfEachLevel(newLevelCallback, getCellContributionWrapper);
  }
   */


  T get_field_chi2() {

    T chi2 = 0;

    for (size_t i = 0; i < getNumLevels(); ++i) {

      auto &field = pGrid[i]->getFieldFourier();
      const auto &spectrum = C0s[i];
      T norm = T(pGrid[i]->size3);

      for (size_t j = 0; j < field.size(); ++j) {
        if (spectrum[j] != 0)
          chi2 += pow(abs(field[j]), 2.0) / (spectrum[j] * norm);
      }
    }

    return chi2;

  }

};

#endif
