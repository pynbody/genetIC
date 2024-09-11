#ifndef _DUMMYIC_HPP_INCLUDED
#define _DUMMYIC_HPP_INCLUDED

template<typename GridDataType>
class ICGenerator;

#include <string>
#include <tuple>
#include <cassert>
#include <functional>
#include <algorithm>
#include <memory>
#include <limits>
#include <iostream>
#include <list>


#include "ic.hpp"
#include "bindings.hpp"



namespace dummyic {

  /*! \class DummyICGenerator
      \brief Special ICGenerator that is used for running parameter files for an input mapper.

      Adds some additional checking when adding levels, to ensure that we use an input mapper
      consistent with the underlying ICGenerator. Also, overrides some functions of the ICgenerator
      class to suppress actions that aren't useful for the purposes of the input mapper.
  */
  template<typename GridDataType>
  class DummyICGenerator : public ICGenerator<GridDataType> {
  protected:
    ICGenerator<GridDataType> *pUnderlying; //!< Underlying ICGenerator calling the input mapper
  public:
    using typename ICGenerator<GridDataType>::T;

    //! Constructor with a specified ICGenerator
    DummyICGenerator(ICGenerator<GridDataType> *pUnderlying) : pUnderlying(pUnderlying) {

    }

    //! Adds a level to this dummy context
    void addLevelToContext(T gridSize, size_t nside, const Coordinate<T> &offset) override {
      if(pUnderlying != nullptr) {
        size_t newLevel = this->multiLevelContext.getNumLevels();
        // getNumLevels counts from 1 to N rather than 0 to N-1, which is why newLevel defined this way does not exist yet
        std::shared_ptr<grids::Grid<T>> underlyingGrid;


        if (pUnderlying->multiLevelContext.getNumLevels() <= newLevel) {
          // source file has extra zoom levels compared to us. Make a grid with our specifications, and any
          // flags deposited onto it will have to be manually copied over later.

          grids::Grid<T> &deepestUnderlyingGrid =
              pUnderlying->multiLevelContext.getGridForLevel(pUnderlying->multiLevelContext.getNumLevels() - 1);

          underlyingGrid = std::make_shared<grids::Grid<T>>(deepestUnderlyingGrid.periodicDomainSize, nside,
                                                            gridSize / nside, offset.x, offset.y, offset.z);
        } else {
          underlyingGrid = pUnderlying->multiLevelContext.getGridForLevel(newLevel).shared_from_this();

        }

        if (underlyingGrid->size != nside)
          throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid n)");

        if (underlyingGrid->thisGridSize != gridSize)
          throw std::runtime_error(
              "Trying to match particles between incompatible simulation setups (wrong grid size)");

        if (!underlyingGrid->offsetLower.almostEqual(offset))
          throw std::runtime_error(
              "Trying to match particles between incompatible simulation setups (wrong grid origin)");

        this->multiLevelContext.addLevel(underlyingGrid);
        this->gadgetTypesForLevels.push_back(1);
      } else {
        ICGenerator<GridDataType>::addLevelToContext(gridSize, nside, offset);
      }
    }

    /* User-level calls from parameter file that require an unnecessary (and potentially costly) operation
     * when working out the relationship between input mappers should not be executed.
     * Ensure this is always the case by overriding them empty in this dummy IC class
    */

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void zeroLevel(size_t /*level*/, size_t) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void dumpGrid(size_t /*level*/, particle::species) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void dumpPS(size_t, particle::species) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void dumpMask() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void dumpID(std::string /*fname*/) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void write() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void modify(std::string /*name*/, std::string /*string*/, float /*value*/) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void done() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void applyModifications() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void calculate(std::string /* name */) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void getFieldChi2() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void reverse() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void reverseSmallK(T /*kmax*/) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void importLevel(size_t /*level*/, std::string /*filename*/) override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void saveTipsyArray(std::string fname, size_t nField) override {}

    /* Override low-levels functions to ensure that certain operations such as drawing the random field or convolving with
     * the power spectrum are never applied in the context of working out relationships between input mappers,
     * covering cases when new user-level facilities are added but not necessarily overriden in the dummy IC class
     * (see https://github.com/pynbody/genetIC/issues/73)
     * */

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void initialiseAllRandomComponents() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void ensureParticleGeneratorInitialised() override {}

    //! Calls to this function has no effect in a dummy IC generator, since it is only working out the mapper structure
    void applyPowerSpec() override {}
  };

  template<typename GridDataType>
  std::shared_ptr<DummyICGenerator<GridDataType>> performDummyRun(const std::string &fname, ICGenerator<GridDataType> *pUnderlying) {
    // Set up the command interpreter to issue commands to main_generator
    tools::ClassDispatch<ICGenerator<GridDataType>, void> dispatch_generator;
    setup_parser(dispatch_generator);

    auto pseudoICs = std::make_shared<DummyICGenerator<GridDataType>>(pUnderlying);

    auto dispatch = dispatch_generator.specify_instance(*pseudoICs);

    std::ifstream inf;
    inf.open(fname);

    if (!inf.is_open())
      throw std::runtime_error("Cannot open IC paramfile for relative_to command");
    logging::entry() << std::endl;
    logging::entry() << "+ Input mapper: computing geometry from " << fname << std::endl;
    {
      tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(fname));
      logging::IndentWhileInScope temporaryIndent;
      logging::entry() << std::endl;
      dispatch.run_loop(inf);
    }
    logging::entry() << std::endl;

    return pseudoICs;

  }

}

#endif
