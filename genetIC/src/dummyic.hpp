#ifndef _DUMMYIC_HPP_INCLUDED
#define _DUMMYIC_HPP_INCLUDED

#include "ic.hpp"

template<typename GridDataType>
class DummyICGenerator : public ICGenerator<GridDataType> {
protected:
  ICGenerator<GridDataType> *pUnderlying;
public:
  using typename ICGenerator<GridDataType>::T;

  DummyICGenerator(ICGenerator<GridDataType> *pUnderlying) : ICGenerator<GridDataType>(pUnderlying->interpreter),
                                                             pUnderlying(pUnderlying) {

  }

  void addLevelToContext(const cosmology::CAMB<GridDataType> & /*&spectrum*/, T gridSize, size_t nside,
                         const Coordinate<T> &offset = {0, 0, 0}) override {
    size_t newLevel = this->multiLevelContext.getNumLevels();
    std::shared_ptr<grids::Grid<T>> underlyingGrid;
    std::shared_ptr<const fields::Field<GridDataType, T>> covarianceField;

    if (pUnderlying->multiLevelContext.getNumLevels() <= newLevel) {
      // source file has extra zoom levels compared to us. Make a grid with our specifications, then
      // make a proxy grid to point into the existing deepest level.
      grids::Grid<T> & deepestUnderlyingGrid =
        pUnderlying->multiLevelContext.getGridForLevel(this->multiLevelContext.getNumLevels()-1);
      grids::Grid<T> gridSpecification(deepestUnderlyingGrid.simsize, nside, gridSize / nside, offset.x, offset.y, offset.z);
      covarianceField = nullptr;
      underlyingGrid = deepestUnderlyingGrid.makeProxyGridToMatch(gridSpecification);
    } else {
      underlyingGrid = pUnderlying->multiLevelContext.getGridForLevel(newLevel).shared_from_this();
      covarianceField = pUnderlying->multiLevelContext.getCovariance(newLevel).shared_from_this();
    }

    if (underlyingGrid->size != nside)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid n)");

    if (underlyingGrid->boxsize != gridSize)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid size)");

    if (underlyingGrid->offsetLower != offset)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid origin)");

    this->multiLevelContext.addLevel(covarianceField,underlyingGrid);
  }


  void zeroLevel(int /*level*/) override {

  }

  void applyPowerSpec() override {}

  void dumpGrid(int /*level*/) override {}

  void dumpPS(int /*level*/) override {}

  virtual void initialiseParticleGenerator() override {}

  void dumpID(string /*fname*/) override {}

  void write() override {}

  void constrain(string /*name*/, string /*string*/, float /*value*/) override {}

  void done() override {}
};

#endif
