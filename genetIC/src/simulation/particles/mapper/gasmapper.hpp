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
      size_t dm0;
      size_t nFirst, nSecond;


      virtual void incrementIterator(iterator *pIterator) const {

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

      virtual size_t size() const {
        return nFirst + nSecond;
      }

      virtual size_t size_gas() const {
        return gasFirst ? nFirst : nSecond;
      }

      virtual size_t size_dm() const {
        return gasFirst ? nSecond : nFirst;
      }

      AddGasMapper(MapPtrType &pFirst, MapPtrType &pSecond, bool gasFirst = true) :
        firstMap(pFirst), secondMap(pSecond), gasFirst(gasFirst), nFirst(pFirst->size()),
        nSecond(pSecond->size()) {
        assert(pFirst->size_gas() == 0);
        assert(pSecond->size_gas() == 0);
      };


      virtual void flagParticles(const std::vector<size_t> &genericParticleArray) {
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
