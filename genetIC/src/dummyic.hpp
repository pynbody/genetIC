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

  void addLevelToContext(const cosmology::CAMB<GridDataType> & spectrum, T gridSize, size_t nside,
                         const Coordinate<T> &offset = {0, 0, 0}) override {
    size_t newLevel = this->multiLevelContext.getNumLevels(); // getNumLevels counts from 1 to N rather than 0 to N-1, which is why newLevel defined this way does not exist yet
    std::shared_ptr<grids::Grid<T>> underlyingGrid;
    std::vector<std::shared_ptr<const fields::Field<GridDataType, T>>> covarianceFieldPtr;


    if (pUnderlying->multiLevelContext.getNumLevels() <= newLevel) {
      // source file has extra zoom levels compared to us. Make a grid with our specifications, and any
      // flags deposited onto it will have to be manually copied over later.

      grids::Grid<T> &deepestUnderlyingGrid =
          pUnderlying->multiLevelContext.getGridForLevel(pUnderlying->multiLevelContext.getNumLevels() - 1);

      underlyingGrid = std::make_shared<grids::Grid<T>>(deepestUnderlyingGrid.periodicDomainSize, nside,
                                                        gridSize / nside, offset.x, offset.y, offset.z);
    } else {
      underlyingGrid = pUnderlying->multiLevelContext.getGridForLevel(newLevel).shared_from_this();

      for(size_t i = 0; i < this->outputFields.size(); i++)
      {
          // TODO - this is very messy - is there a better way?
        try {
        auto resPointer = pUnderlying->multiLevelContext.getCovariance(newLevel,i).shared_from_this();
        covarianceFieldPtr.push_back(resPointer);
        } catch (const std::out_of_range &e) {
              // leave covarianceFieldPtr as nullptr
              covarianceFieldPtr.push_back(nullptr);
        }
      }
    }

    if (underlyingGrid->size != nside)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid n)");

    if (underlyingGrid->thisGridSize != gridSize)
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid size)");

    if (!underlyingGrid->offsetLower.almostEqual(offset))
      throw std::runtime_error("Trying to match particles between incompatible simulation setups (wrong grid origin)");

    this->multiLevelContext.addLevel(covarianceFieldPtr, underlyingGrid);
    this->gadgetTypesForLevels.push_back(1);
  }


  void zeroLevel(size_t /*level*/,size_t) override {}

  void applyPowerSpec(size_t) override {}

  void dumpGrid(size_t /*level*/,size_t) override {}

  void dumpPS(size_t,size_t) override {}

  void dumpMask() override {}

  virtual void initialiseParticleGenerator(size_t) override {}

  void dumpID(string /*fname*/) override {}

  void write() override {}

  void modify(string /*name*/, string /*string*/, float /*value*/) override {}

  void done() override {}
};

#endif
