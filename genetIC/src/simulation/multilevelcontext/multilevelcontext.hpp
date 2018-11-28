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
    std::vector<std::vector<std::shared_ptr<const fields::Field<DataType, T>>>> C0s;
    std::vector<T> weights;

  public:
    size_t nTransferFunctions;


  protected:
    std::vector<size_t> Ns;
    std::vector<size_t> cumu_Ns;
    size_t Ntot;
    size_t nLevels;
    T simSize;

    MultiLevelContextInformationBase(size_t N) {
      nLevels = 1;
      Ns.push_back(N);
      Ntot = N;
      nTransferFunctions = 1;
    }


  public:

    MultiLevelContextInformationBase() {
      nLevels = 0;
      Ntot = 0;
      nTransferFunctions = 1;
    }

    virtual ~MultiLevelContextInformationBase() {}


    void
    addLevel(const cosmology::CAMB<DataType> &spectrum, T size, size_t nside, const Coordinate<T> &offset = {0, 0, 0}) {

      if (!spectrum.isUsable())
        throw std::runtime_error("Cannot add a grid level until the power spectrum has been specified");

      if (nLevels == 0)
        simSize = size;

      auto grid = std::make_shared<grids::Grid<T>>(simSize, nside, size / nside, offset.x, offset.y, offset.z);

      nTransferFunctions = spectrum.dmOnly ? 1 : spectrum.nTransfers;

      std::vector<std::shared_ptr<const fields::Field<DataType, T>>> C0 (nTransferFunctions);
      for(size_t i = 0;i < nTransferFunctions;i++)
      {
        C0[i] = spectrum.getPowerSpectrumForGrid(*grid,i);
      }
      addLevel(C0, grid);
    }

    void addLevel(const std::vector<std::shared_ptr<const fields::Field<DataType, T>>>& C0, std::shared_ptr<grids::Grid<T>> pG) {
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

    // For adding a vector where each element is the same:
    void addLevel(std::shared_ptr<const fields::Field<DataType, T>> C0,std::shared_ptr<grids::Grid<T>> pG)
    {
        std::vector<std::shared_ptr<const fields::Field<DataType, T>>> constArray (this->nTransferFunctions,C0);
        this->addLevel(constArray,pG);
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

    //! \brief Return a reference to the grid vector. Needed to add gas to all levels, for example.
    std::vector<std::shared_ptr<grids::Grid<T>>> getAllGrids()
    {
        return pGrid;
    }

    T getWeightForLevel(size_t level) {
      return weights[level];
    }

    // Main function used by the code to retrieve the power spectrum.
    // TODO - find a better way than this cascade of ifs...
    const fields::Field<DataType> &getCovariance(size_t level,size_t nTransfer = 0) const {
      if (C0s.size() > level)
      {
          if(C0s[level].size() > nTransfer)
          {
              if(C0s[level][nTransfer] != nullptr)
                {
                    return *C0s[level][nTransfer];
                }
          }
      }
      // Throw if not returned anything:
      throw std::out_of_range("No covariance yet specified for this level");
    }

    size_t deepestLevelwithFlaggedCells() {

      for (int i = this->getNumLevels() - 1; i >= 0; --i) {
        if (this->getGridForLevel(i).numFlaggedCells()>0)
          return size_t(i);
      }

      throw std::runtime_error("No level has any particles selected");
    }

    size_t deepestCoarseGrid() const {
      for (int i = this->getNumLevels() - 1; i >= 0; --i) {
        if (this->getGridForLevel(i).coversFullSimulation())
          return size_t(i);
      }
      throw std::runtime_error("No coarse grid were found.");
    }

    //! From finest level, use interpolation to construct other levels
    //! Since multiLevelContext doesn't know which type of field it is creating
    //! (baryons or dark matter) we have to supply this information as a parameter.
    //! \param data - field data
    //! \param nField - nField = 0 for dark matter (default), nField = 1 for baryons.
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
                                                    size_t extra_lores, size_t extra_highres) const {
      //! Copy this MultiLevelContextInformation, but insert intermediate virtual grids such that
      //!there is a full stack increasing in the specified power.
      /*!
       * E.g. if there is a 256^3 and a 1024^3 grid stack, with default parameters this will return a
       * 128^3, 256^3, 512^3 and 1024^3 grid stack.
       *
       * The first two will be based on the 256^3 literal grid, and the second two on the 1024^3 literal grid.
       *
       * @param newStack     the MultiLevelContextInformation into which the new stack will be placed. Any
       *                     existing grids in the stack will be removed.
       * @param base_factor  grids will be downgraded by factors of base_factor^N where N is an integer
       * @param extra_lores  number of additional grids *above* the base level to add
       * @param extra_highes  number of additional grids *below* the finest level to add
      */
      newStack.clear();

      // TODO Refactor this loop to make it neater.
      // First intermediate resolution
      // Second Low res subsample stuffed from the base level
      // Third High res supersampled stuff from the finest level
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

      size_t factor = base_factor;
      for (size_t i = 0; i < extra_highres; ++i) {
        size_t level = nLevels - 1;
        size_t neff = size_t(round(pGrid[0]->cellSize / pGrid[level]->cellSize)) * pGrid[0]->size;
        std::cerr << "Adding virtual grid with effective resolution " << neff * factor << std::endl;
        auto vGrid = std::make_shared<grids::SuperSampleGrid<T>>(pGrid[level], factor);
        newStack.addLevel(nullptr, vGrid);
        factor *= base_factor;

      }

    }

    void copyContextAndCenter(MultiLevelContextInformation<DataType> &newStack,
                              const Coordinate<T> pointToCenterOnto) const {

      newStack.clear();
      size_t deepest_coarse = this->deepestCoarseGrid();
      assert(deepest_coarse < nLevels);
      Coordinate<T> offset;

      for (size_t level = 0; level <= deepest_coarse; ++level) {
        auto centeredCoarse = std::make_shared<grids::CenteredGrid<T>>(this->pGrid[level], pointToCenterOnto);
        newStack.addLevel(C0s[level], centeredCoarse);
        offset = centeredCoarse->getPointOffset();
      }


      for (size_t level = deepest_coarse + 1; level < nLevels; ++level) {
        auto offsetFine = std::make_shared<grids::OffsetGrid<T>>(this->pGrid[level], offset.x, offset.y, offset.z);
        newStack.addLevel(C0s[level], offsetFine);
      }
    }

    void copyContextWithCenteredIntermediate(MultiLevelContextInformation<DataType> &newStack,
                                             const Coordinate<T> pointToCenterOnto,
                                             size_t base_factor,
                                             size_t subsample, size_t supersample) const {

      // Transform resolution ratios into number of grids
      int extralowres;
      int extrahighres;
      try {
        extralowres = tools::findPowerOf(base_factor, subsample);
        extrahighres = tools::findPowerOf(base_factor, supersample);
      } catch (const std::runtime_error &e){
        throw std::runtime_error("Subsample and supersample cannot be converted into an integer number of grids with this base factor");
      }

      auto extracontext = multilevelcontext::MultiLevelContextInformation<DataType>();
      this->copyContextWithIntermediateResolutionGrids(extracontext, base_factor, extralowres, extrahighres);
      extracontext.copyContextAndCenter(newStack, pointToCenterOnto);
    }


    size_t getIndexOfCellOnOtherLevel(size_t currentLevel, size_t otherLevel, size_t cellIndex){
      Coordinate<T> cell_coord(this->getGridForLevel(currentLevel).getCentroidFromIndex(cellIndex));
      return this->getGridForLevel(otherLevel).getIndexFromPoint(cell_coord);
    }

    std::vector<size_t> getFullResolutionGrids(){
      std::vector<size_t> levelsOfRealGrids;

      for (size_t i = 0; i<nLevels; i++){
        if(! this->getGridForLevel(i).isUpsampledOrDownsampled()){
          levelsOfRealGrids.push_back(i);
        }
      }
      assert(! levelsOfRealGrids.empty());
      return levelsOfRealGrids;
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
};

#endif
