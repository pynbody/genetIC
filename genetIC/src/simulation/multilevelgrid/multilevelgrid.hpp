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
#include "src/simulation/particles/species.hpp"

namespace fields {
  template<typename T>
  class ConstraintField;
}

namespace cosmology {
  template<typename T>
  class PowerSpectrum;
}

namespace filters {
  template<typename T>
  class FilterFamily;
}

/*!
    \namespace multilevelgrid
    \brief Keep track of the different grids and their level of refinement. Intermediate between grids and multi level fields.

    Because C++ doesn't allow partial specialization of member functions, there is a workaround involving splitting
    most of the class into a Base class, and the final specializations are made in a child class.
    See http://stackoverflow.com/questions/165101/invalid-use-of-incomplete-type-error-with-partial-template-specialization

 */
namespace multilevelgrid {

  /*! \class MultiLevelGrid
  \brief Object used to store information about the grids on each zoom level

  Inherits from the base class in order to specialise to real/complex grids.
*/
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class MultiLevelGrid;

  /*! \class MultiLevelGridBase
    \brief Base class defining the multi-level grid. Contains most of the class structure
  */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class MultiLevelGridBase : public tools::Signaling {
  private:
    std::vector<std::shared_ptr<grids::Grid<T>>> pGrids; //!< Pointers to the grids for each level
    std::vector<std::shared_ptr<grids::Grid<T>>> pOutputGrid; //!< Pointers to the output grids for each level -- may be different from the underlying grids if allowStrays is on
    std::vector<T> weights; //!< Fraction of the volume of the coarsest level's cells that the cells on each level occupy
    const cosmology::PowerSpectrum<DataType> *powerSpectrumGenerator = nullptr;

    bool levelsAreCombined = false; //!< Set to true after low-frequency information is propagated into high-resolution fields

  public:
    size_t nTransferFunctions = 1; //!< Keeps track of the number of transfer functions currently being used by the code.
    bool allowStrays = false; //!< If true, return output grids that cover the whole simulation box even in the zoom regions

  protected:
    std::vector<size_t> Ns; //!< Vector that stores the number of cells on each level
    std::vector<size_t> cumu_Ns; //!< Vector that stores the cumulative number of cells, ie, for each level, holds the number of cells on that level and all levels above it
    size_t Ntot = 0; //!< Total number of cells in the multi-level context
    size_t nLevels = 0; //!< Number of levels in the multi-level context.
    T simSize; //!< Comoving size of the simulation box


    MultiLevelGridBase() {}


  public:


    virtual ~MultiLevelGridBase() {}

    void setLevelsAreCombined() {
      levelsAreCombined = true;
    }

    bool getLevelsAreCombined() const {
      return levelsAreCombined;
    }


    /*! \brief Adds a level to the current multi-level context
        \param size - co-moving size of the level to be added
        \param nside - number of cells on one side of the grid of the level being added
        \param offset - offset of the level being added from the lower front left hand corner of the coarsest grid. Defaults to (0,0,0)
    */
    void
    addLevel(T size, size_t nside, const Coordinate<T> &offset = {0, 0, 0}) {

      if (nLevels == 0)
        simSize = size;

      auto grid = std::make_shared<grids::Grid<T>>(simSize, nside, size / nside, offset.x, offset.y, offset.z);

      nTransferFunctions = 2;
      addLevel(grid);
    }

    /*! Define the power spectrum generator to be used when power spectra are required.
     *
     * This must remain alive for the lifetime of the MultiLevelGrid class
    */
    void setPowerspectrumGenerator(const cosmology::PowerSpectrum<DataType> &generator) {
      this->powerSpectrumGenerator = &generator;
    }

    /*! \brief Performs the process of actually adding the level with the appropriate grid and transfer functions
        \param pG - pointer to the grid for this level
    */
    void addLevel(std::shared_ptr<grids::Grid<T>> pG) {
      if (pGrids.size() == 0) {
        weights.push_back(1.0);
      } else {
        weights.push_back(pow(pG->cellSize / pGrids[0]->cellSize, 3.0));
      }
      pGrids.push_back(pG);
      Ns.push_back(pG->size3);
      cumu_Ns.push_back(Ntot);
      Ntot += pG->size3;

      if (allowStrays && nLevels > 0) {
        // output grid covers entire simulation box, using interpolation from next level up to fill the gaps
        auto gridForOutput = std::make_shared<grids::ResolutionMatchingGrid<T>>(pG, pOutputGrid[nLevels - 1]);
        pOutputGrid.push_back(gridForOutput);
      } else {
        // output grid is identical to grid on which data is actually stored
        pOutputGrid.push_back(pG);
      }
      nLevels += 1;
      this->changed();
    }

    //! Clears the multi-level context of all data
    void clear() {
      pGrids.clear();
      Ns.clear();
      cumu_Ns.clear();
      Ntot = 0;
      nLevels = 0;
      levelsAreCombined = false;
      this->changed();
    }


    //! Returns the total number of cells in the multi-level context
    size_t getNumCells() const {
      return Ntot;
    }

    //! Returns the number of degrees of freedom -- i.e. the number of cells, counting only the finest in any location
    size_t getNumDof() const {
      size_t dof = 0;
      for(size_t i=0; i<getNumLevels(); ++i) {
        dof+=getGridForLevel(i).size3;
        if(i>0) {
          // subtract the number of cells this corresponds to on the level above
          size_t factor = tools::getRatioAndAssertInteger(getGridForLevel(i-1).cellSize,getGridForLevel(i).cellSize);
          factor*=factor*factor;
          dof-=getGridForLevel(i).size3/factor;
        }
      }
      return dof;
    }

    //! Returns the number of levels in the multi-level context
    size_t getNumLevels() const {
      return nLevels;
    }

    //! Returns a reference to the calculation grid (on which data is manipulated) for the specified level
    grids::Grid<T> &getGridForLevel(size_t level) {
      return *pGrids[level];
    }

    //! Returns a constant reference to the calculation grid (on which data is manipulated) for the specified level
    const grids::Grid<T> &getGridForLevel(size_t level) const {
      return *pGrids[level];
    }

    //! Returns a reference to the output grid (from which particles are output) for the specified level
    grids::Grid<T> &getOutputGridForLevel(size_t level) {
      return *pOutputGrid[level];
    }

    //! Returns a constant reference to the output grid (from which particles are output) for the specified level
    const grids::Grid<T> &getOutputGridForLevel(size_t level) const {
      return *pOutputGrid[level];
    }

    //! Return a vector of shared pointers to the grids. Needed to add gas to all levels, for example.
    std::vector<std::shared_ptr<grids::Grid<T>>> getAllGrids() {
      return pGrids;
    }

    //! Returns the fractional volume of a level-0 pixel that is occupied by a pixel on this level
    T getWeightForLevel(size_t level) const {
      return weights[level];
    }

    /*! \brief Returns a reference to the specified transfer function for the specified level
        \param level - level to get transfer function for
        \param species - the type of particle, which will potentially determine which transfer function is used
    */
    std::shared_ptr<fields::Field<DataType>> getCovariance(size_t level, particle::species species) const {
      // Caching is now implemented in the power spectrum, so copies of the power spectrum on each level are no
      // longer stored in this class.
      assert(this->powerSpectrumGenerator);
      return this->powerSpectrumGenerator->getPowerSpectrumForGrid(this->getGridForLevel(level).shared_from_this(),
                                                                   species);
    }

    //! Returns the index of the deepest level that has flagged cells
    size_t deepestLevelwithFlaggedCells() {

      for (int i = this->getNumLevels() - 1; i >= 0; --i) {
        if (this->getGridForLevel(i).numFlaggedCells() > 0)
          return size_t(i);
      }

      throw std::runtime_error("No level has any particles selected");
    }

    //! Returns the index of the deepest 'coarse grid', ie, the deepest level that covers the entire simulation.
    size_t deepestCoarseGrid() const {
      for (int i = this->getNumLevels() - 1; i >= 0; --i) {
        if (this->getGridForLevel(i).coversFullSimulation())
          return size_t(i);
      }
      throw std::runtime_error("No coarse grid were found.");
    }

    /*! \brief From finest level, use interpolation to construct other levels

        \param data - field data on the highest resolution level.
    */
    std::shared_ptr<fields::ConstraintField<DataType>>
    generateMultilevelCovectorFromHiresCovector(fields::Field<DataType, T> &&data, particle::species transferType) const {
      using tools::numerics::operator*=;
      assert(&data.getGrid() == pGrids.back().get());

      // Generate the fields on each level. Fill low-res levels with zeros to start with.
      vector<std::shared_ptr<fields::Field<DataType, T>>> fieldsOnLevels;
      for (size_t level = 0; level < pGrids.size(); level++) {
        if (level == pGrids.size() - 1) {
          fieldsOnLevels.emplace_back(std::make_shared<fields::Field<DataType, T>>(std::move(data)));
        } else {
          fieldsOnLevels.emplace_back(std::make_shared<fields::Field<DataType, T>>(*pGrids[level], false));
        }
      }

      // Now interpolate the high-res level down into lower-res levels
      size_t levelmax = fieldsOnLevels.size() - 1;
      if (levelmax > 0) {
        assert(fieldsOnLevels.back()->isFourier());
        fieldsOnLevels.back()->toReal();

        for (int level = levelmax - 1; level >= 0; --level) {

          fields::Field<DataType, T> & hires = *fieldsOnLevels[level + 1];
          fields::Field<DataType, T> & lores = *fieldsOnLevels[level];

          int pixel_size_ratio = tools::getRatioAndAssertInteger(fieldsOnLevels[level]->getGrid().cellSize,
                                                                 fieldsOnLevels[level+1]->getGrid().cellSize);
          int pixel_volume_ratio = pow(pixel_size_ratio,3);



          /* Naive loop that can't be parallelized:
          size_t hiresFieldSize = hires.getGrid().size3;
          for(size_t i=0; i<hiresFieldSize; i++) {
            auto location = hires.getGrid().getCentroidFromIndex(i);
            lores.deInterpolate(location, hires[i]/pixel_volume_ratio);
          }
          */

          // deinterpolation affects a grid of 4^3 cells (on the lores grid), so to parallelise we need to ensure
          // that threads are always working on regions at least 4 lores cells apart...
          size_t hiresGridSize = hires.getGrid().size;
          size_t hiresParallelisationSeparation = 4 * pixel_size_ratio;
          for(int xoffset=0; xoffset < hiresParallelisationSeparation; xoffset++) {
#pragma omp parallel for default(none) shared(hires, lores, pixel_volume_ratio, xoffset, hiresGridSize, hiresParallelisationSeparation)
            for(int x=xoffset; x<hiresGridSize; x+=hiresParallelisationSeparation) {
              for(int y=0; y<hiresGridSize; y++) {
                for(int z=0; z<hiresGridSize; z++) {
                  Coordinate<int> coord {x,y,z};
                  auto index = hires.getGrid().getIndexFromCoordinateNoWrap(coord);
                  auto location = hires.getGrid().getCentroidFromCoordinate(coord);
                  lores.deInterpolate(location, hires[index]/pixel_volume_ratio);
                }
              }
            }
          }





        }
      }

      for(size_t level=0; level < pGrids.size(); ++level) {
        fieldsOnLevels[level]->getDataVector()*= getWeightForLevel(level) / getWeightForLevel(pGrids.size() - 1);
      }

      auto retVal = std::make_shared<fields::ConstraintField<DataType>>(
        *dynamic_cast<const MultiLevelGrid<DataType, T> *>(this),
        fieldsOnLevels, transferType, true);

      retVal->applyFilters(filters::FilterFamily<DataType>(*this));
      retVal->toReal();
      return retVal;
    }

    //! Applies the specified operation to the grids on each level of the multi-level context
    void forEachLevel(std::function<void(grids::Grid<T> &)> newLevelCallback) {
      for (size_t level = 0; level < nLevels; level++) {
        newLevelCallback(*pGrids[level]);
      }
    }

    //! Applies the specified operation (that cannot change anything) to the grids for each level.
    void forEachLevel(std::function<void(const grids::Grid<T> &)> newLevelCallback) const {
      for (size_t level = 0; level < nLevels; level++) {
        newLevelCallback(*pGrids[level]);
      }
    }


    void copyContextWithIntermediateResolutionGrids(MultiLevelGrid<DataType> &newStack,
                                                    size_t resolution_step_factor = 2,
                                                    size_t number_extra_lowres = 0,
                                                    size_t number_extra_highres = 0) const {
      //! Copy this MultiLevelGrid, but insert intermediate virtual grids such that
      //! there is a full stack increasing in the specified power.
      /*!
      * E.g. if there is a 256^3 and a 1024^3 grid stack, with default parameters this will return a
      * 128^3, 256^3, 512^3 and 1024^3 grid stack.
      *
      * The first two will be based on the 256^3 literal grid, and the second two on the 1024^3 literal grid.
      *
      * @param newStack     the MultiLevelGrid into which the new stack will be placed. Any
      *                     existing grids in the stack will be removed.
      * @param resolution_step_factor   grids will be downgraded by factors of base_factor^N where N is an integer
      * @param number_extra_lowres  number of additional grids *above* the base level to add
      * @param number_extra_highres  number of additional grids *below* the finest level to add
      */

      newStack.clear();

      for (size_t level = 0; level < nLevels; ++level) {
        if (level == 0) {
          // First add the base subsampled grids from 'subsample' command
          this->addExtraSubSampledGrids(newStack, number_extra_lowres, resolution_step_factor, level);

          logging::entry() << "Adding real grid with resolution " << pGrids[0]->size << std::endl;
          newStack.addLevel(pGrids[level]);
        } else if (level > 0) {

          // Add the intermediate grids between base and zoom grids
          size_t neff_zoom = size_t(round(pGrids[0]->cellSize / pGrids[level]->cellSize)) * pGrids[0]->size;
          double zoomfactor = (pGrids[level - 1]->cellSize / pGrids[level]->cellSize);

          auto number_intermediate_grids = size_t(
            round(log10(zoomfactor) / log10(double(resolution_step_factor))) - 1);
          this->addExtraSubSampledGrids(newStack, number_intermediate_grids, resolution_step_factor, level);

          logging::entry() << "Adding real grid with resolution " << neff_zoom << std::endl;
          newStack.addLevel(pGrids[level]);
        }
      }
      // Add supersampled grids from the 'supersample' command
      this->addExtraSuperSampledGrids(newStack, number_extra_highres, resolution_step_factor);
    }

    void addExtraSubSampledGrids(MultiLevelGridBase<DataType> &newStack,
                                 size_t number_extra_grids, size_t resolution_step_factor, size_t level) const {


      auto factor = size_t(pow(resolution_step_factor, number_extra_grids));
      size_t neff_real_grid = size_t(round(pGrids[0]->cellSize / pGrids[level]->cellSize)) * pGrids[0]->size;

      if (number_extra_grids > 0) {
        for (size_t i = 0; i < number_extra_grids; ++i) {
          auto vGrid = std::make_shared<grids::SubSampleGrid<T>>(pGrids[level], factor);

          logging::entry() << "Adding virtual grid with effective resolution " << neff_real_grid / factor << std::endl;
          newStack.addLevel(vGrid);

          factor /= resolution_step_factor;
        }
      }
    }

    void addExtraSuperSampledGrids(MultiLevelGridBase<DataType> &newStack,
                                   size_t number_extra_highres, size_t resolution_step_factor) const {

      size_t factor = resolution_step_factor;
      size_t finest_level = nLevels - 1;
      size_t neff_finest = size_t(round(pGrids[0]->cellSize / pGrids[finest_level]->cellSize)) * pGrids[0]->size;

      if (number_extra_highres > 0) {
        for (size_t i = 0; i < number_extra_highres; ++i) {
          auto vGrid = std::make_shared<grids::SuperSampleGrid<T>>(pGrids[finest_level], factor);

          logging::entry() << "Adding virtual grid with effective resolution " << neff_finest * factor << std::endl;
          newStack.addLevel(vGrid);

          factor *= resolution_step_factor;
        }
      }
    }

    /*! \brief Creates a copy of the multi-level context, but centred on the specified point.
        \param newStack - reference to where the new context should be stored
        \param pointToCenterOnto - new centre point.
    */
    void copyContextAndCenter(MultiLevelGrid<DataType> &newStack,
                              const Coordinate<T> pointToCenterOnto) const {

      newStack.clear();
      size_t deepest_coarse = this->deepestCoarseGrid();
      assert(deepest_coarse < nLevels);
      Coordinate<T> offset;

      for (size_t level = 0; level <= deepest_coarse; ++level) {
        auto centeredCoarse = std::make_shared<grids::CenteredGrid<T>>(this->pGrids[level], pointToCenterOnto);
        newStack.addLevel(centeredCoarse);
        offset = centeredCoarse->getPointOffset();
      }


      for (size_t level = deepest_coarse + 1; level < nLevels; ++level) {
        auto offsetFine = std::make_shared<grids::OffsetGrid<T>>(this->pGrids[level], offset.x, offset.y, offset.z);
        newStack.addLevel(offsetFine);
      }
    }

    /*! \brief Copies the multi-level context, adding intermediate resolution grids, and centering them all on the specified point.
        \param newStack - reference to where the new context should be stored
        \param pointToCenterOnto - new centre point.
        \param base_factor - grids will be downgraded by factors of base_factor^N where N is an integer
        \param subsample - factor times smaller resolution that we want the coarsest grid to be
        \param supersample - factor times higher resolution that we want the finest grid to be
    */
    void copyContextWithCenteredIntermediate(MultiLevelGrid<DataType> &newStack,
                                             const Coordinate<T> pointToCenterOnto,
                                             size_t base_factor,
                                             size_t subsample, size_t supersample) const {

      // Transform resolution ratios into number of grids
      int extralowres;
      int extrahighres;
      try {
        extralowres = tools::findPowerOf(base_factor, subsample);
        extrahighres = tools::findPowerOf(base_factor, supersample);
      } catch (const std::runtime_error &e) {
        throw std::runtime_error(
          "Subsample and supersample cannot be converted into an integer number of grids with this base factor");
      }

      auto extracontext = multilevelgrid::MultiLevelGrid<DataType>();
      this->copyContextWithIntermediateResolutionGrids(extracontext, base_factor, extralowres, extrahighres);
      extracontext.copyContextAndCenter(newStack, pointToCenterOnto);
    }


    /*! \brief Returns the index that the cell on the specified level corresponds to on another level
        \param currentLevel - level the specified index corresponds to currently
        \param otherLevel - level we want the index for corresponding to the input index
        \param cellIndex - index on the current level
    */
    size_t getIndexOfCellOnOtherLevel(size_t currentLevel, size_t otherLevel, size_t cellIndex) {
      Coordinate<T> cell_coord(this->getGridForLevel(currentLevel).getCentroidFromIndex(cellIndex));
      return this->getGridForLevel(otherLevel).getIndexFromPoint(cell_coord);
    }

    //! Returns a vector of all the levels that are neither upsampled nor downsampled.
    std::vector<size_t> getFullResolutionGrids() {
      std::vector<size_t> levelsOfRealGrids;

      for (size_t i = 0; i < nLevels; i++) {
        if (!this->getGridForLevel(i).isUpsampledOrDownsampled()) {
          levelsOfRealGrids.push_back(i);
        }
      }
      assert(!levelsOfRealGrids.empty());
      return levelsOfRealGrids;
    }

  };

  /*! \class MultiLevelGrid<T, T>
    \brief Object used to store information about the multi-level context

    Inherits from the base class in order to specialise to real/complex grids.
  */
  template<typename T>
  class MultiLevelGrid<T, T> : public MultiLevelGridBase<T, T> {
  protected:
    using DataType = T;
    using MultiLevelGridBase<DataType, T>::cumu_Ns;
    using MultiLevelGridBase<DataType, T>::Ns;
    using MultiLevelGridBase<DataType, T>::nLevels;

  public:
    /*! \brief Sums the result of evaluating the specified function over the cells on the specified levels
        \param newLevelCallback - vector of bools indicating which levels to include in the sum
        \param getCellContribution - function to evaluate on each cell
    */
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
