#ifndef IC_ONELEVELMAPPER_HPP
#define IC_ONELEVELMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"

namespace particle {

  namespace mapper {
    using std::endl;
    using std::cerr;


    template<typename GridDataType>
    class OneLevelParticleMapper : public ParticleMapper<GridDataType> {

    private:

      using MapType = ParticleMapper<GridDataType>;
      using typename MapType::T;
      using typename MapType::MapPtrType;
      using typename MapType::GridPtrType;
      using typename MapType::ConstGridPtrType;
      using typename MapType::GridType;
      using typename MapType::iterator;

      GridPtrType pGrid;

      unsigned int gadgetParticleType;

    protected:

      virtual void
      dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
        gp = pGrid;
        i = pIterator->i;
      }

      virtual void decrementIteratorBy(iterator *pIterator, size_t increment) const override {
        pIterator->i -= increment;
      }

      virtual unsigned int
      gadgetParticleTypeFromIterator(const iterator * /*pIterator*/) const override {
        return gadgetParticleType;
      }

    public:
      virtual void debugInfo(std::ostream &s, int level = 0) const override {
        tools::indent(s, level);
        pGrid->debugInfo(s);
        s << std::endl;
      }

      void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const override {
        tools::indent(s, n);
        s << "i=" << pIterator->i << " into iterator for grid " << pGrid << endl;
      }

      OneLevelParticleMapper(std::shared_ptr<grids::Grid<T>> &pGrid) : pGrid(pGrid) {
        gadgetParticleType = 1;
      }

      OneLevelParticleMapper(std::shared_ptr<grids::Grid<T>> &&pGrid) : pGrid(pGrid) {
        gadgetParticleType = 1;
      }

      size_t size() const override {
        return pGrid->size3;
      }

      virtual void unflagAllParticles() override {
        pGrid->unflagAllCells();
      }

      void flagParticles(const std::vector<size_t> &genericParticleArray) override {
        pGrid->flagCells(genericParticleArray);
      }

      void getFlaggedParticles(std::vector<size_t> &particleArray) const override {
        pGrid->getFlaggedCells(particleArray);
      }

      GridPtrType getCoarsestGrid() override {
        return pGrid;
      }

      GridPtrType getFinestGrid() override {
        return pGrid;
      }

      bool supportsReverseIterator() override {
        return true;
      }

      void setGadgetParticleType(unsigned int type) {
        gadgetParticleType = type;
      }

      virtual iterator beginParticleType(const AbstractMultiLevelParticleGenerator<GridDataType> &generator,
                                         unsigned int particleType) const override {
        if(gadgetParticleType==particleType) {
          return this->begin(generator);
        } else {
          return iterator(nullptr, generator);
        }
      }

      virtual iterator endParticleType(const AbstractMultiLevelParticleGenerator<GridDataType> &generator,
                                       unsigned int particleType) const override {
        if(gadgetParticleType==particleType) {
          return this->end(generator);
        } else {
          return iterator(nullptr, generator);
        }
      }

      std::pair<MapPtrType, MapPtrType> addGas(T massRatio, const std::vector<GridPtrType> &toGrids) override {
        if (pGrid->isProxyForAnyOf(toGrids)) {
          return std::make_pair(
              std::make_shared<OneLevelParticleMapper<GridDataType>>(
                  this->pGrid->makeScaledMassVersion(massRatio)),
              std::make_shared<OneLevelParticleMapper<GridDataType>>(
                  this->pGrid->makeScaledMassVersion(1.0 - massRatio)));
        } else {
          return std::make_pair(std::shared_ptr<OneLevelParticleMapper<GridDataType>>(nullptr),
                                std::make_shared<OneLevelParticleMapper<GridDataType>>(this->pGrid));
        }
      }

      MapPtrType superOrSubSampleDM(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {
        std::shared_ptr<OneLevelParticleMapper<GridDataType>> newMapper;
        if (pGrid->isProxyForAnyOf(toGrids)) {
          GridPtrType newGrid;
          if (super)
            newGrid = std::make_shared<grids::SuperSampleGrid<T >>(this->pGrid, ratio);
          else
            newGrid = std::make_shared<grids::SubSampleGrid<T >>(this->pGrid, ratio);
          newMapper = std::make_shared<OneLevelParticleMapper<GridDataType>>(newGrid);
        } else {
          newMapper = std::make_shared<OneLevelParticleMapper<GridDataType>>(this->pGrid);
        }
        newMapper->setGadgetParticleType(this->gadgetParticleType);
        return newMapper;
      }


    };
  }
}
#endif
