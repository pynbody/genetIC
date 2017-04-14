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

  void addLevelToContext(const cosmology::CAMB<T> & /*&spectrum*/, T gridSize, size_t nside,
                         const Coordinate<T> &offset = {0, 0, 0}) override {
    size_t newLevel = this->multiLevelContext.getNumLevels();
    if (pUnderlying->multiLevelContext.getNumLevels() <= newLevel)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (too many levels)");

    grids::Grid<T> &underlyingGrid = pUnderlying->multiLevelContext.getGridForLevel(newLevel);
    if (underlyingGrid.size != nside)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid n)");

    if (underlyingGrid.boxsize != gridSize)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid size)");

    if (underlyingGrid.offsetLower != offset)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid origin)");

    this->multiLevelContext.addLevel(pUnderlying->multiLevelContext.getCovariance(newLevel),
                                     underlyingGrid.shared_from_this());
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
