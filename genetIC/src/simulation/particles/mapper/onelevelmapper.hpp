#ifndef IC_ONELEVELMAPPER_HPP
#define IC_ONELEVELMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"

namespace particle {

  namespace mapper {
    using std::endl;
    using std::cerr;


    /*! \class OneLevelParticleMapper
        \brief Mapper used to refer to the particles on a single level of the multi-level context.

        This is the base object used to construct two-level mappers. Each one level mapper deals with
        particles on a single cubic grid, without any knowledge of other levels that may exist.
    */
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

      GridPtrType pGrid; //!< Pointer to the grid corresponding to the level for this mapper

      unsigned int gadgetParticleType; //!< Gadget type for the particle species pointed to by this mapper

    protected:

      //! Dereferences the iterator with the grid for this level, and the current position
      virtual void
      dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
        gp = pGrid;
        i = pIterator->i;
      }

      /*! \brief Decrements the iterator by the specified step.
        \param pIterator - iterator to decrement.
        \param increment - number of steps to decrement the iterator by.
      */
      virtual void decrementIteratorBy(iterator *pIterator, size_t increment) const override {
        pIterator->i -= increment;
      }

      //! Returns the gadget type associated to this one level mapper
      virtual unsigned int
      gadgetParticleTypeFromIterator(const iterator * /*pIterator*/) const override {
        return gadgetParticleType;
      }

    public:
      /*! \brief Outputs debug information about the one level mapper to the specified stream
        \param s - stream to output to.
        \param level - level about which to output debug information.
      */
      virtual void debugInfo(std::ostream &s, int level = 0) const override {
        tools::indent(s, level);
        pGrid->debugInfo(s);
        s << std::endl;
      }

      /*! \brief Outputs debug information about the specified iterator
      */
      void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const override {
        tools::indent(s, n);
        s << "i=" << pIterator->i << " into iterator for grid " << pGrid << endl;
      }

      //! Constructor with a single grid corresponding to the level this is defined on.
      OneLevelParticleMapper(std::shared_ptr<grids::Grid<T>> &pGrid) : pGrid(pGrid) {
        gadgetParticleType = 1;
      }

      //! Constructor with a single grid corresponding to the level this is defined on.
      OneLevelParticleMapper(std::shared_ptr<grids::Grid<T>> &&pGrid) : pGrid(pGrid) {
        gadgetParticleType = 1;
      }

      //! Returns the number of cells at this level
      size_t size() const override {
        return pGrid->size3;
      }

      //! Unflags all the flagged cells on this level
      virtual void unflagAllParticles() override {
        pGrid->unflagAllCells();
      }

      //! Flags the specified cells on this grid
      void flagParticles(const std::vector<size_t> &genericParticleArray) override {
        pGrid->flagCells(genericParticleArray);
      }

      //! Copies the flagged cells to the specified vector
      void getFlaggedParticles(std::vector<size_t> &particleArray) const override {
        pGrid->getFlaggedCells(particleArray);
        std::sort(particleArray.begin(), particleArray.end());
      }

      //! Returns the grid associated to this level (only one level, so by definition it is the coarsest)
      GridPtrType getCoarsestGrid() override {
        return pGrid;
      }

      //! Returns the grid associated to this level (only one level, so by definition it is the finest)
      GridPtrType getFinestGrid() override {
        return pGrid;
      }

      //! One level mappers support reverse iteration
      bool supportsReverseIterator() override {
        return true;
      }

      //! Change the gadget type of this one level mapper (default is dark matter)
      void setGadgetParticleType(unsigned int type) {
        gadgetParticleType = type;
      }

      //! Go to the begining of the particles of the specified gadget type.
      virtual iterator beginParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                         unsigned int particleType) const override {
        if (gadgetParticleType == particleType) {
          return this->begin(generator);
        } else {
          return iterator(nullptr, generator);
        }
      }

      //! Go to the end of the particles of the specified gadget type.
      virtual iterator endParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                       unsigned int particleType) const override {
        if (gadgetParticleType == particleType) {
          return this->end(generator);
        } else {
          return iterator(nullptr, generator);
        }
      }

      /*! \brief Creates a gas mapper from this one level mapper, out of two mass scaled grids with the appropriate baryon and dark matter mass fractions.
        \param massRatio - fraction of the mass in a given particle species
        \param toGrids - vector of grids - only create usable mappers if this grid is a proxy for one of these.
      */
      std::pair<MapPtrType, MapPtrType> splitMass(T massRatio, const std::vector<GridPtrType> &toGrids) override {
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

      /*! \brief sub/super samples dark matter
        \param ratio - factor to sub/super sample by
        \param toGrids - sub/super-sample if at least one of these grids is a proxy for this one
        \param super - if true, super-samples the dark matter, otherwise the dark matter is sub-sampled.
      */
      MapPtrType superOrSubSample(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {
        std::shared_ptr<OneLevelParticleMapper<GridDataType>> newMapper;
        if (pGrid->isProxyForAnyOf(toGrids)) {
          GridPtrType newGrid;
          if (super)
            newGrid = this->pGrid->makeSupersampled(ratio);
          else
            newGrid = this->pGrid->makeSubsampled(ratio);
          newMapper = std::make_shared<OneLevelParticleMapper<GridDataType>>(newGrid);
        } else {
          newMapper = std::make_shared<OneLevelParticleMapper<GridDataType>>(this->pGrid);
        }
        newMapper->setGadgetParticleType(this->gadgetParticleType);
        return newMapper;
      }

      MapPtrType insertIntermediateResolutionPadding(size_t, size_t) override {
        return this->shared_from_this();
      }

      MapPtrType withIndependentFlags() override {
        return std::make_shared<OneLevelParticleMapper<GridDataType>>(this->pGrid->withIndependentFlags());
      }

      MapPtrType withCoupledFlags() override {
        return std::make_shared<OneLevelParticleMapper<GridDataType>>(this->pGrid->withCoupledFlags());
      }


    };
  }
}
#endif
