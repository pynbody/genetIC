#ifndef IC_GASMAPPER_HPP
#define IC_GASMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"

namespace particle {

  namespace mapper {
    using std::endl;
    using std::cerr;


    /*! \class AddGasMapper
        \brief Ties together a pair of mappers, one for gas (baryons), one for dark matter.

        Essentially, this object consists of a pair of mappers, one which points to dark matter,
        and another which points to baryons. Which is which depends on the order in which
        particles are required in the output format. Operations applied to particles are then
        propagated down to the appropriate mapper.
    */
    template<typename GridDataType>
    class AddGasMapper : public ParticleMapper<GridDataType> {
    public:
      using MapType = ParticleMapper<GridDataType>;
      using typename MapType::T;
      using typename MapType::MapPtrType;
      using typename MapType::GridPtrType;
      using typename MapType::GridType;
      using typename MapType::iterator;
      using typename MapType::ConstGridPtrType;

    protected:
      MapPtrType firstMap; //!< First mapper in the pair. Whether this is baryons or dark matter depends on the output format
      MapPtrType secondMap; //!< Second mapper in the pair
      bool gasFirst; //!< If true, baryons are assigned using the first mapper in the pair, and if false, using the second.
      size_t nFirst; //!< Number of particles in the first mapper.
      size_t nSecond; //!< Number of particles in the second mapper.


      //! Increments the specified iterator, increasing either its baryon or dark matter iterator according to the current position
      virtual void incrementIterator(iterator *pIterator) const override {

        if (pIterator->i >= nFirst)
          ++(*(pIterator->subIterators[1]));
        else
          ++(*(pIterator->subIterators[0]));
        ++(pIterator->i);

      }


      /*! \brief Dereferences either the first or second iterator according to the current position, returning the grid and particle index pointed to
        \param pIterator - constant pointer to the iterator to de-reference
        \param gp - reference to where the resulting grid pointer should be stored
        \param i - reference to where the resulting cell index should be stored
      */
      virtual void
      dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
        if (pIterator->i >= nFirst)
          pIterator->subIterators[1]->deReference(gp, i);
        else
          pIterator->subIterators[0]->deReference(gp, i);
      }


    public:


      //! Returns true if either the first or second mapper points to the specified grid
      bool references(GridPtrType grid) const override {
        return firstMap->references(grid) || secondMap->references(grid);
      }

      /*! \brief Outputs debug information about the specified level to the specified stream.
        \param s - stream to output debug information to.
        \param level - level to output debug information about.
      */
      virtual void debugInfo(std::ostream &s, int level = 0) const override {
        tools::indent(s, level);
        s << "AddGasMapper";
        if (gasFirst)
          s << ", gas first:" << std::endl;
        else
          s << ", DM first:" << std::endl;

        firstMap->debugInfo(s, level + 1);

        tools::indent(s, level);

        s << "AddGasMapper";
        if (gasFirst)
          s << ", DM second:" << std::endl;
        else
          s << ", gas second:" << std::endl;

        secondMap->debugInfo(s, level + 1);

        tools::indent(s, level);
        s << "AddGasMapper ends" << std::endl;
      }

      //! Returns the total number of particles in the two mappers
      virtual size_t size() const override {
        return nFirst + nSecond;
      }

      //! Returns the number of baryons on the baryon mappers
      virtual size_t size_gas() const override {
        return gasFirst ? nFirst : nSecond;
      }

      //! Returns the number of dark matter particles in the dark matter mapper
      virtual size_t size_dm() const override {
        return gasFirst ? nSecond : nFirst;
      }

      /*! \brief Constructor taking the two mappers to be bound, and the order to assign them in
        \param pFirst - pointer to first mapper
        \param pSecond - pointer to the second mapper
        \param gasFirst - if true, pFirst is baryons, pSecond is dark matter, and vice versa if false
      */
      AddGasMapper(MapPtrType pFirst, MapPtrType pSecond, bool gasFirst = true) :
        firstMap(pFirst), secondMap(pSecond), gasFirst(gasFirst), nFirst(pFirst->size()),
        nSecond(pSecond->size()) {
        assert(pFirst->size_gas() == 0);
        assert(pSecond->size_gas() == 0);
        for (unsigned int particleType = 0; particleType < 6; ++particleType) {
          // All particle types must be the default (1), as we will now override them for our own purposes
          if (particleType != 1) {
            NullMultiLevelParticleGenerator<GridDataType> g;
            if (pFirst->beginParticleType(g, particleType) != pFirst->endParticleType(g, particleType) ||
                pSecond->beginParticleType(g, particleType) != pSecond->endParticleType(g, particleType)) {
              throw std::runtime_error("Cannot currently combine custom gadget particle type numbers with gas");
            }
          }
        }
      };


      //! Flag particles specified by genericParticleArray (currently only flags dark matter particles)
      virtual void flagParticles(const std::vector<size_t> &genericParticleArray) override {
        std::vector<size_t> first;
        std::vector<size_t> second;

        for (auto i: genericParticleArray) {
          if (i < nFirst)
            first.push_back(i);
          else
            second.push_back(i - nFirst);
        }

        // At present we ONLY distribute particles onto the DM grids
        if (gasFirst)
          secondMap->flagParticles(second);
        else
          firstMap->flagParticles(first);

      }

      //! Unflags all the particles in both mappers
      virtual void unflagAllParticles() override {
        firstMap->unflagAllParticles();
        secondMap->unflagAllParticles();
      }

      //! Returns the coarsest grid associated to the dark matter mapper
      virtual GridPtrType getCoarsestGrid() override {
        return gasFirst ? secondMap->getCoarsestGrid() : firstMap->getCoarsestGrid();
      }

      //! Returns the finest grid associated to the dark matter mapper
      GridPtrType getFinestGrid() override {
        return gasFirst ? secondMap->getFinestGrid() : firstMap->getFinestGrid();
      }

      //! Copies the flags (currently only dark matter flags) from the underlying grids to the particleArray vector
      void getFlaggedParticles(std::vector<size_t> &particleArray) const override {
        // At present we ONLY gather particles from the DM grids, to make
        // this the inverse operation to flagParticles.
        if (gasFirst) {
          secondMap->getFlaggedParticles(particleArray);

          // perform offset
          for (size_t &i : particleArray)
            i += nFirst;

        } else {
          firstMap->getFlaggedParticles(particleArray);
        }
      }

      /*! \brief Return an iterator pointing to the beginning of the list of particles of the specified type
        \param generator - particle generator to use
        \param particleType - particle species to find the beginning of (0 for dark matter, 1 for baryons)
      */
      virtual iterator beginParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                         unsigned int particleType) const override {
        if (particleType == 0)
          return beginGas(generator);
        else if (particleType == 1)
          return beginDm(generator);
        else
          return endDm(generator);
      }


      /*! \brief Return an iterator pointing to the end of the list of particles of the specified type
        \param generator - particle generator to use
        \param particleType - particle species to find the end of (0 for dark matter, 1 for baryons)
      */
      virtual iterator endParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                       unsigned int particleType) const override {
        if (particleType == 0)
          return endGas(generator);
        else if (particleType == 1)
          return endDm(generator);
        else
          return endDm(generator);
      }


      //! Returns an iterator set to the beginning of the particle list
      virtual iterator begin(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        iterator i(this, generator);
        i.subIterators.emplace_back(new iterator(firstMap->begin(generator)));
        i.subIterators.emplace_back(new iterator(secondMap->begin(generator)));
        return i;
      }

      //! Returns an iterator pointing to the beginning of the dark matter particles
      virtual iterator
      beginDm(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? secondMap : firstMap)->begin(generator);
      }

      //! Returns an iterator pointing to the end fo the dark matter particles
      virtual iterator endDm(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? secondMap : firstMap)->end(generator);
      }

      //! Returns an iterator pointing to the beginning of the baryonic particles
      virtual iterator
      beginGas(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? firstMap : secondMap)->begin(generator);
      }

      //! Returns an iterator pointing to the end of the dark matter particles
      virtual iterator endGas(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? firstMap : secondMap)->end(generator);
      }

      /*! \brief Returns a mapper that either super-samples or sub-samples the dark matter
        \param ratio - ratio for sub/super-samples to use
        \param toGrids - vector of pointers to the grids to sub/super sample.
        \param super - if true, super-samples the dark matter. Otherwise, subsamples it.
      */
      MapPtrType superOrSubSample(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {
        auto ssub1 = firstMap;
        auto ssub2 = secondMap;

        bool applyTo2 = gasFirst || (!super);
        bool applyTo1 = (!gasFirst) || (!super);

        if (applyTo2)
          ssub2 = ssub2->superOrSubSample(ratio, toGrids, super);
        if (applyTo1)
          ssub1 = ssub1->superOrSubSample(ratio, toGrids, super);

        return std::make_shared<AddGasMapper<GridDataType>>(
          ssub1, ssub2, gasFirst);
      }

      MapPtrType insertIntermediateResolutionPadding(size_t ratio, size_t padCells) override {
        auto ssub1 = firstMap;
        auto ssub2 = secondMap;

        if (gasFirst)
          ssub2 = ssub2->insertIntermediateResolutionPadding(ratio, padCells);
        else
          ssub1 = ssub1->insertIntermediateResolutionPadding(ratio, padCells);

        return std::make_shared<AddGasMapper<GridDataType>>(
          ssub1, ssub2, gasFirst);
      }

      MapPtrType withIndependentFlags() override {
        return std::make_shared<AddGasMapper<T>>(firstMap->withIndependentFlags(), secondMap->withIndependentFlags(),
                                                 gasFirst);
      }

      MapPtrType withCoupledFlags() override {
        return std::make_shared<AddGasMapper<T>>(firstMap->withCoupledFlags(), secondMap->withCoupledFlags(), gasFirst);
      }

    };
  }
}
#endif
