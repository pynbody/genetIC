#ifndef IC_GASMAPPER_HPP
#define IC_GASMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"

namespace particle {

  namespace mapper {
    using std::endl;
    using std::cerr;


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
      MapPtrType firstMap;
      MapPtrType secondMap;
      bool gasFirst;
      size_t nFirst, nSecond;


      virtual void incrementIterator(iterator *pIterator) const override {

        if (pIterator->i >= nFirst)
          ++(*(pIterator->subIterators[1]));
        else
          ++(*(pIterator->subIterators[0]));
        ++(pIterator->i);

      }


      virtual void
      dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
        if (pIterator->i >= nFirst)
          pIterator->subIterators[1]->deReference(gp, i);
        else
          pIterator->subIterators[0]->deReference(gp, i);
      }


    public:



      bool references(GridPtrType grid) const override {
        return firstMap->references(grid) || secondMap->references(grid);
      }

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

      virtual size_t size() const override {
        return nFirst + nSecond;
      }

      virtual size_t size_gas() const override {
        return gasFirst ? nFirst : nSecond;
      }

      virtual size_t size_dm() const override {
        return gasFirst ? nSecond : nFirst;
      }

      AddGasMapper(MapPtrType &pFirst, MapPtrType &pSecond, bool gasFirst = true) :
          firstMap(pFirst), secondMap(pSecond), gasFirst(gasFirst), nFirst(pFirst->size()),
          nSecond(pSecond->size()) {
        assert(pFirst->size_gas() == 0);
        assert(pSecond->size_gas() == 0);
        for(unsigned int particleType=0; particleType<6; ++particleType) {
          // All particle types must be the default (1), as we will now override them for our own purposes
          if(particleType!=1) {
            NullMultiLevelParticleGenerator<GridDataType> g;
            if(pFirst->beginParticleType(g,particleType)!=pFirst->endParticleType(g,particleType) ||
              pSecond->beginParticleType(g,particleType)!=pSecond->endParticleType(g,particleType)) {
              throw std::runtime_error("Cannot currently combine custom gadget particle type numbers with gas");
            }
          }
        }
      };


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
        // MR: Not sure what this means?
        if (first.empty() && !second.empty())
          secondMap->flagParticles(second);
        else
          firstMap->flagParticles(first);

      }

      virtual void unflagAllParticles() override {
        firstMap->unflagAllParticles();
        secondMap->unflagAllParticles();
      }

      virtual GridPtrType getCoarsestGrid() override {
        return gasFirst ? secondMap->getCoarsestGrid() : firstMap->getCoarsestGrid();
      }

      GridPtrType getFinestGrid() override {
        return gasFirst ? secondMap->getFinestGrid() : firstMap->getFinestGrid();
      }

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

      virtual iterator beginParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                         unsigned int particleType) const override {
        if(particleType==0)
          return beginGas(generator);
        else if(particleType==1)
          return beginDm(generator);
        else
          return endDm(generator);
      }


      virtual iterator endParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                       unsigned int particleType) const override {
        if(particleType==0)
          return endGas(generator);
        else if(particleType==1)
          return endDm(generator);
        else
          return endDm(generator);
      }


      virtual iterator begin(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        iterator i(this, generator);
        i.subIterators.emplace_back(new iterator(firstMap->begin(generator)));
        i.subIterators.emplace_back(new iterator(secondMap->begin(generator)));
        return i;
      }

      virtual iterator
      beginDm(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? secondMap : firstMap)->begin(generator);
      }

      virtual iterator endDm(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? secondMap : firstMap)->end(generator);
      }

      virtual iterator
      beginGas(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? firstMap : secondMap)->begin(generator);
      }

      virtual iterator endGas(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return (gasFirst ? firstMap : secondMap)->end(generator);
      }

      MapPtrType superOrSubSampleDM(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {
        auto ssub1 = firstMap;
        auto ssub2 = secondMap;

        bool applyTo2 = gasFirst || (!super);
        bool applyTo1 = (!gasFirst) || (!super);

        if (applyTo2)
          ssub2 = ssub2->superOrSubSampleDM(ratio, toGrids, super);
        if (applyTo1)
          ssub1 = ssub1->superOrSubSampleDM(ratio, toGrids, super);

        return std::make_shared<AddGasMapper<GridDataType>>(
            ssub1, ssub2, gasFirst);
      }

    };
  }
}
#endif
