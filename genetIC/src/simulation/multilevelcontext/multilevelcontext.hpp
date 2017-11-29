#ifndef _MULTILEVELFIELDMANAGER_HPP
#define _MULTILEVELFIELDMANAGER_HPP

#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <cassert>
#include <memory>

#include "src/simulation/grid/grid.hpp"
#include "src/tools/signaling.hpp"


namespace fields {
  template<typename T>
  class ConstraintField;
}

namespace cosmology {
  template<typename T>
  class CAMB;
}

/*!
    \namespace multilevelcontext
    \brief Keep track of the different grids and their level of refinement. Intermediate between grids and multi level fields.

    Because C++ doesn't allow partial specialization of member functions, there is a workaround involving splitting
    most of the class into a Base class, and the final specializations are made in a child class.
    See http://stackoverflow.com/questions/165101/invalid-use-of-incomplete-type-error-with-partial-template-specialization

 */
namespace multilevelcontext {
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class MultiLevelContextInformation;

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class MultiLevelContextInformationBase : public tools::Signaling {
  private:
    std::vector<std::shared_ptr<grids::Grid<T>>> pGrid;
    std::vector<std::shared_ptr<const fields::Field<DataType, T>>> C0s;
    std::vector<T> weights;

  protected:
    std::vector<size_t> Ns;
    std::vector<size_t> cumu_Ns;
    size_t Ntot;
    size_t nLevels;
    T simSize;

    //TODO Never used, should delete ?
//    void mapIdToLevelId(size_t i, size_t &level, size_t &level_id) {
//      level = 0;
//      while (level < nLevels && i >= Ns[level]) {
//        i -= Ns[level];
//        level++;
//      }
//      if (i >= Ns[level])
//        throw std::runtime_error("ID out of range when mapping into underlying fields in mapIdToLevelId");
//      level_id = i;
//    }

    MultiLevelContextInformationBase(size_t N) {
      nLevels = 1;
      Ns.push_back(N);
      Ntot = N;
    }


  public:

    MultiLevelContextInformationBase() {
      nLevels = 0;
      Ntot = 0;
    }

    virtual ~MultiLevelContextInformationBase() {}


    void
    addLevel(const cosmology::CAMB<DataType> &spectrum, T size, size_t nside, const Coordinate<T> &offset = {0, 0, 0}) {
      if (!spectrum.isUsable())
        throw std::runtime_error("Cannot add a grid level until the power spectrum has been specified");

      if (nLevels == 0)
        simSize = size;

      auto grid = std::make_shared<grids::Grid<T>>(simSize, nside, size / nside, offset.x, offset.y, offset.z);
      auto C0 = spectrum.getPowerSpectrumForGrid(*grid);
      addLevel(C0, grid);

    }

    void addLevel(std::shared_ptr<const fields::Field<DataType, T>> C0, std::shared_ptr<grids::Grid<T>> pG) {
      C0s.push_back(C0);
      if (pGrid.size() == 0) {
        weights.push_back(1.0);
      } else {
        weights.push_back(pow(pG->cellSize / pGrid[0]->cellSize, 3.0));
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

    grids::Grid<T> &getGridForLevel(size_t level) {
      return *pGrid[level];
    }

    const grids::Grid<T> &getGridForLevel(size_t level) const {
      return *pGrid[level];
    }

    T getWeightForLevel(size_t level) {
      return weights[level];
    }

    const fields::Field<DataType> &getCovariance(size_t level) const {
      if (C0s.size() > level && C0s[level] != nullptr)
        return *C0s[level];
      else
        throw std::out_of_range("No covariance yet specified for this level");
    }

    size_t deepestLevelwithFlaggedCells() {

      for (int i = this->getNumLevels() - 1; i >= 0; --i) {
        if (this->getGridForLevel(i).hasFlaggedCells())
          return size_t(i);
      }

      throw std::runtime_error("No level has any particles selected");
    }

    //! From finest level, use interpolation to construct other levels
    std::shared_ptr<fields::ConstraintField<DataType>>
    generateMultilevelFromHighResField(fields::Field<DataType, T> &&data) {
      assert(&data.getGrid() == pGrid.back().get());

      // Generate the fields on each level. Fill low-res levels with zeros to start with.
      vector<std::shared_ptr<fields::Field<DataType, T>>> dataOnLevels;
      for (size_t level = 0; level < pGrid.size(); level++) {
        if (level == pGrid.size() - 1) {
          dataOnLevels.emplace_back(std::make_shared<fields::Field<DataType, T>>(std::move(data)));
        } else {
          dataOnLevels.emplace_back(std::make_shared<fields::Field<DataType, T>>(*pGrid[level], false));
        }
      }

      // Now interpolate the high-res level down into lower-res levels
      size_t levelmax = dataOnLevels.size() - 1;
      if (levelmax > 0) {
        assert(dataOnLevels.back()->isFourier());
        dataOnLevels.back()->toReal();
        for (int level = levelmax - 1; level >= 0; --level) {
          dataOnLevels[level]->addFieldFromDifferentGrid(*(dataOnLevels[level + 1]));
        }
      }

      return std::make_shared<fields::ConstraintField<DataType>>(
          *dynamic_cast<MultiLevelContextInformation<DataType, T> *>(this),
          dataOnLevels);
    }

    void forEachLevel(std::function<void(grids::Grid<T> &)> newLevelCallback) {
      for (size_t level = 0; level < nLevels; level++) {
        newLevelCallback(*pGrid[level]);
      }
    }

    void forEachLevel(std::function<void(const grids::Grid<T> &)> newLevelCallback) const {
      for (size_t level = 0; level < nLevels; level++) {
        newLevelCallback(*pGrid[level]);
      }
    }


    void copyContextWithIntermediateResolutionGrids(MultiLevelContextInformation<DataType> &newStack,
                                                    size_t base_factor,
                                                    size_t extra_lores) const {
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
        size_t neff = size_t(round(pGrid[0]->cellSize / pGrid[level]->cellSize)) * pGrid[0]->size;
        if (level > 0) {
          size_t factor = base_factor;
          while (pGrid[level]->cellSize * factor * 1.001 < pGrid[level - 1]->cellSize) {
            std::cerr << "Adding virtual grid with effective resolution " << neff / factor << std::endl;
            auto vGrid = std::make_shared<grids::SubSampleGrid<T>>(pGrid[level], factor);
            newStack.addLevel(nullptr, vGrid);
            factor *= base_factor;
          }
        } else {
          size_t factor = base_factor;
          for (size_t i = 0; i < extra_lores; ++i) {
            std::cerr << "Adding virtual grid with effective resolution " << neff / factor << std::endl;
            auto vGrid = std::make_shared<grids::SubSampleGrid<T>>(pGrid[level], factor);
            newStack.addLevel(nullptr, vGrid);
            factor *= base_factor;
          }
        }

        std::cerr << "Adding real grid with resolution " << neff << std::endl;
        newStack.addLevel(C0s[level], pGrid[level]);

      }

    }

    size_t getIndexOfCellOnOtherLevel(size_t currentLevel, size_t otherLevel, size_t cellIndex){
      auto currentLevelGrid = this->getGridForLevel(currentLevel);
      auto otherLevelGrid = this->getGridForLevel(otherLevel);
      Coordinate<T> cell_coord(currentLevelGrid.getCellCentroid(cellIndex));
      return otherLevelGrid.getCellContainingPoint(cell_coord);
    }

  };

  template<typename T>
  class MultiLevelContextInformation<T, T> : public MultiLevelContextInformationBase<T, T> {
  protected:
    using DataType = T;
    using MultiLevelContextInformationBase<DataType, T>::cumu_Ns;
    using MultiLevelContextInformationBase<DataType, T>::Ns;
    using MultiLevelContextInformationBase<DataType, T>::nLevels;

  public:
    DataType accumulateOverEachModeOfEachLevel(
        std::function<bool(size_t)> newLevelCallback,
        std::function<DataType(size_t, size_t, size_t)> getCellContribution) {

      // real version only - see complex specialization below

      T accum(0);

      for (size_t level = 0; level < nLevels; level++) {

        if (!newLevelCallback(level))
          continue;

        size_t level_base = cumu_Ns[level];

#pragma omp parallel for reduction(+:accum)
        for (size_t i = 0; i < Ns[level]; i++) {
          size_t i_all_levels = level_base + i;
          DataType res = getCellContribution(level, i, i_all_levels);

          // accumulate separately - OMP doesn't support complex number reduction :-(
          accum += res;
        }
      }

      return accum;
    }
  };

  //TODO Is this class used at all ? It looks an unused specialization of base class
  template<typename T>
  class MultiLevelContextInformation<std::complex<T>, T>
      : public MultiLevelContextInformationBase<std::complex<T>, T> {
  protected:
    using DataType= std::complex<T>;
    using MultiLevelContextInformationBase<DataType, T>::cumu_Ns;
    using MultiLevelContextInformationBase<DataType, T>::Ns;
    using MultiLevelContextInformationBase<DataType, T>::nLevels;

  public:


  };
};

#endif
