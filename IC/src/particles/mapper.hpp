//
// mapper.hpp
//
// The idea of the classes in this file are to provide a way to map particle IDs
// onto grid locations. This becomes quite complicated when one has multiple
// grids and potentially things like gas particles and so on. By making
// specific classes to do very specific bits of the mapping, hopefully we keep
// this complexity managable while also having a lot of flexibility to handle
// different set-ups.

#ifndef __MAPPER_HPP
#define __MAPPER_HPP

#include <memory>
#include <typeinfo>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "src/grid.hpp"
#include "src/argsort.hpp"
#include "src/multilevelfield.hpp"
#include "multilevelgenerator.hpp"

// helper function for our debug dumps
void indent(std::ostream &s, int level = 0) {
  for (int i = 0; i < level; i++) s << "| ";
  s << "+ ";

}

namespace particle {

// forward declarations:
  template<typename T>
  class ParticleMapper;

  template<typename T>
  class OneLevelParticleMapper;

  template<typename T>
  class TwoLevelParticleMapper;

  template<typename T>
  class AddGasMapper;

  template<typename T>
  class ParticleGenerator;

  template<typename T>
  class AbstractMultiLevelParticleGenerator;

  template<typename T>
  class MapperIterator {
  protected:
    // This is a generic structure used by all subclasses to keep
    // track of where they are in a sequential operation

    using MapType = ParticleMapper<T>;
    using MapPtrType = std::shared_ptr<ParticleMapper<T>>;
    using GridType = Grid<T>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;
    using DereferenceType = std::pair<ConstGridPtrType, size_t>;

    friend class ParticleMapper<T>;

    friend class OneLevelParticleMapper<T>;

    friend class TwoLevelParticleMapper<T>;

    friend class AddGasMapper<T>;



    size_t i;
    std::vector<std::shared_ptr<MapperIterator<T>>> subIterators;
    std::vector<size_t> extraData;
    const ParticleMapper<T> *pMapper;
    const AbstractMultiLevelParticleGenerator<T> &generator;
    mutable ConstGridPtrType pLastGrid;
    mutable std::shared_ptr<const ParticleGenerator<T>> pLastGridGenerator;

    MapperIterator(const ParticleMapper<T> *pMapper,
                   const AbstractMultiLevelParticleGenerator<T> &generator) :
            i(0), pMapper(pMapper), generator(generator) {}


  public:

    MapperIterator(const MapperIterator<T> &source) :
          i(source.i), extraData(source.extraData), pMapper(source.pMapper),
          pLastGrid(nullptr), generator(source.generator)
    {
      for (const auto &subIterator: source.subIterators) {
        if (subIterator == nullptr)
          subIterators.push_back(nullptr);
        else
          subIterators.push_back(std::make_shared<MapperIterator<T>>(*subIterator));
      }
    }

    MapperIterator &operator++() {
      pMapper->incrementIterator(this);
      return (*this);
    }

    MapperIterator operator++(int) {
      MapperIterator R = (*this);
      pMapper->incrementIterator(this);
      return R;
    }

    MapperIterator &operator+=(size_t m) {
      pMapper->incrementIteratorBy(this, m);
      return (*this);
    }

    MapperIterator &operator-=(size_t m) {
      pMapper->decrementIteratorBy(this, m);
      return (*this);
    }

    DereferenceType operator*() const {
      ConstGridPtrType gp;
      size_t i;
      deReference(gp, i);
      return std::make_pair(gp, i);
    }


    void deReference(ConstGridPtrType &gp, size_t &id) const {
      // TODO: weird that this is now partly updating internal state and partly returning results
      pMapper->dereferenceIterator(this, gp, id);
      if(gp!=pLastGrid) {
        pLastGrid = gp;
        updateGridReference();
      }
    }


    particle::Particle<T> getParticle() const {
      ConstGridPtrType pGrid;
      size_t id;
      deReference(pGrid, id);
      return pLastGridGenerator->getParticle(*pGrid, id);
    }

  protected:
    void updateGridReference() const {
      pLastGridGenerator = generator.getGeneratorForGrid(*pLastGrid).shared_from_this();
    }

  public:

    template<typename S>
    auto getField(const MultiLevelField<S> &multiLevelField) const {
      const auto q = **this;
      return multiLevelField.getFieldForGrid(*q.first)[q.second];
    }


    size_t getNextNParticles(std::vector<particle::Particle<T>> &particles) {
      size_t n = 1024 * 256;
      if (n + i > pMapper->size())
        n = pMapper->size() - i;

      particles.resize(n);

      n = parallelIterate([&](size_t local_i, const MapperIterator &localIterator) {
        particles[local_i] = localIterator.getParticle();

      }, n);

      return n;
    }

  protected:
    size_t parallelIterate(std::function<void(size_t, const MapperIterator &)> callback, size_t nMax) {
      size_t n = std::min(pMapper->size() - i, nMax);
      size_t final_i = i + n;


#pragma omp parallel
      {
        MapperIterator *pThreadLocalIterator;
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
#else
        int thread_num = 0;
        int num_threads = 1;
#endif

        if (thread_num == 0)
          pThreadLocalIterator = this;
        else
          pThreadLocalIterator = new MapperIterator(*this);

#pragma omp barrier
        (*pThreadLocalIterator) += thread_num;


        for (size_t local_i = thread_num; local_i < n; local_i += num_threads) {
          callback(local_i, *pThreadLocalIterator);
          if (local_i + num_threads < n) (*pThreadLocalIterator) += num_threads;
        }

        if (thread_num != 0)
          delete pThreadLocalIterator;

      }
      if (final_i > i)
        (*this) += final_i - i;
      return n;
    }


  public:

    T getMass() {
      ConstGridPtrType pGrid;
      size_t id;
      deReference(pGrid, id);
      return pLastGridGenerator->getMass(*pGrid);
    }

    std::unique_ptr<DereferenceType> operator->() const {
      return std::unique_ptr<DereferenceType>(new DereferenceType(**this));
    }

    friend bool operator==(const MapperIterator<T> &lhs, const MapperIterator<T> &rhs) {
      return (lhs.i == rhs.i) && (lhs.pMapper == rhs.pMapper);
    }

    friend bool operator!=(const MapperIterator<T> &lhs, const MapperIterator<T> &rhs) {
      return (lhs.i != rhs.i) || (lhs.pMapper != rhs.pMapper);
    }

    void debugInfo() const {
      debugInfo(cerr, 0);
    }

    void debugInfo(std::ostream &s, int n = 0) const {
      pMapper->debugInfoForIterator(s, n, this);
    }


  };


  template<typename T>
  class ParticleMapper {
  public:
    using MapType = ParticleMapper<T>;
    using MapPtrType = std::shared_ptr<ParticleMapper<T>>;
    using GridType = Grid<T>;
    using GridPtrType = std::shared_ptr<Grid<T>>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;
    using iterator = MapperIterator<T>;
    using BaseGeneratorType = particle::AbstractMultiLevelParticleGenerator<T>;

    friend class MapperIterator<T>;

  protected:

    std::shared_ptr<BaseGeneratorType> pGenerator;

    virtual void incrementIterator(iterator *pIterator) const {
      ++(pIterator->i);
    }

    virtual void incrementIteratorBy(iterator *pIterator, size_t increment) const {
      pIterator->i += increment;
    }

    virtual void decrementIteratorBy(iterator *pIterator, size_t increment) const {
      throw std::runtime_error("Attempting to reverse in a mapper that does not support random access");
    }

    virtual void dereferenceIterator(const iterator *pIterator, ConstGridPtrType &, size_t &) const {
      throw std::runtime_error("There is no grid associated with this particle mapper");
    }

  public:

    virtual void debugInfo(std::ostream &s, int level = 0) const {
      indent(s, level);
      s << "Abstract MapperIterator (shouldn't be here)" << std::endl;
    }

    virtual void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const {
      indent(s, n);
      s << "i=" << pIterator->i << " into abstract MapperIterator (\?\?)" << endl;
      for (auto q: pIterator->extraData) {
        indent(s, n);
        s << "data=" << q << endl;
      }
      for (auto q: pIterator->subIterators) {
        indent(s, n);
        if (q != nullptr) {
          s << "subiterator: " << endl;
          q->debugInfo(s, n + 1);
        } else {
          s << "null subiterator" << endl;
        }
      }
    }

    friend std::ostream &operator<<(std::ostream &stream, const ParticleMapper<T> &I) {
      I.debugInfo(stream);
      return stream;
    }


    virtual size_t size() const {
      return 0;
    }

    virtual size_t size_gas() const {
      return 0;
    }

    virtual size_t size_dm() const {
      return size();
    }

    virtual void clearParticleList() {

    }

    virtual void flagParticles(const std::vector<size_t> &genericParticleArray) {
      throw std::runtime_error("Cannot interpret particles yet; no particle->grid mapper available");
    }

    virtual void flagParticles(std::vector<size_t> &&genericParticleArray) {
      // Sometimes it's helpful to have move semantics, but in general just call
      // the normal implementation
      this->flagParticles(genericParticleArray);
    }

    virtual void getFlaggedParticles(std::vector<size_t> &particleArray) const {
      throw std::runtime_error("Cannot get particles; no particle->grid mapper available");
    }

    virtual void extendParticleListToUnreferencedGrids(const std::vector<GridPtrType> &grids) {
      /* For any grid that is _not_ referenced by this mapper, generate cell flags by matching
       * to the finest level available in this mapper.
       *
       * This is used when specifying constraints wrt an unzoomed simulation in a
       * zoomed simulation - the cell flags will not reach onto the finest level until this routine
       * is run.
       */
      for (auto pGrid: grids) {
        if (!this->references(pGrid)) {
          vector<size_t> ar;
          GridPtrType proxyGrid = getFinestGrid()->makeProxyGridToMatch(*pGrid);
          proxyGrid->getFlaggedCells(ar);
          pGrid->flagCells(ar);
        }
      }
    }

    virtual GridPtrType getCoarsestGrid() {
      throw std::runtime_error("There is no grid associated with this particle mapper");
    }

    ConstGridPtrType getCoarsestGrid() const {
      return (const_cast<ParticleMapper *>(this)->getCoarsestGrid());
    }


    virtual GridPtrType getFinestGrid() {
      throw std::runtime_error("There is no grid associated with this particle mapper");
    }

    ConstGridPtrType getFinestGrid() const {
      return (const_cast<ParticleMapper *>(this)->getFinestGrid());
    }

    virtual iterator begin(const AbstractMultiLevelParticleGenerator<T> &generator) const {
      return iterator(this, generator);
    }

    virtual bool references(GridPtrType grid) const {
      return getCoarsestGrid().get() == grid.get();
    }

    virtual iterator end(const AbstractMultiLevelParticleGenerator<T> &generator) const {
      iterator x(this, generator);
      x.i = size();

      return x;
    }

    virtual iterator beginDm(const AbstractMultiLevelParticleGenerator<T> &generator) const {
      return begin(generator);
    }

    virtual iterator endDm(const AbstractMultiLevelParticleGenerator<T> &generator) const {
      return end(generator);
    }

    virtual iterator beginGas(const AbstractMultiLevelParticleGenerator<T> &generator) const {
      return end(generator);
    }

    virtual iterator endGas(const AbstractMultiLevelParticleGenerator<T> &generator) const {
      return end(generator);
    }


    virtual bool supportsReverseIterator() {
      return false;
    }

    virtual std::pair<MapPtrType, MapPtrType> addGas(T massRatio, const std::vector<GridPtrType> &toGrids) {
      throw std::runtime_error("Don't know how to add gas in this context");
    }

    virtual MapPtrType superOrSubSampleDM(int ratio, const std::vector<GridPtrType> &toGrids, bool super = true) {
      throw std::runtime_error("Don't know how to supersample in this context");
    }


  };


  template<typename T>
  class OneLevelParticleMapper : public ParticleMapper<T> {

  private:

    using MapType = ParticleMapper<T>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::ConstGridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;

    GridPtrType pGrid;

  protected:

    virtual void dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
      gp = pGrid;
      i = pIterator->i;
    }

    virtual void decrementIteratorBy(iterator *pIterator, size_t increment) const override {
      pIterator->i -= increment;
    }


  public:
    virtual void debugInfo(std::ostream &s, int level = 0) const override {
      indent(s, level);
      pGrid->debugInfo(s);
      s << std::endl;
    }

    void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const override {
      indent(s, n);
      s << "i=" << pIterator->i << " into iterator for grid " << pGrid << endl;
    }

    OneLevelParticleMapper(std::shared_ptr<Grid<T>> &pGrid) : pGrid(pGrid) {

    }

    OneLevelParticleMapper(std::shared_ptr<Grid<T>> &&pGrid) : pGrid(pGrid) {

    }

    size_t size() const override {
      return pGrid->size3;
    }

    virtual void clearParticleList() override {
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

    std::pair<MapPtrType, MapPtrType> addGas(T massRatio, const std::vector<GridPtrType> &toGrids) override {
      if (pGrid->pointsToAnyGrid(toGrids)) {
        return std::make_pair(
          std::make_shared<OneLevelParticleMapper<T>>(this->pGrid->makeScaledMassVersion(massRatio)),
          std::make_shared<OneLevelParticleMapper<T>>(
            this->pGrid->makeScaledMassVersion(1.0 - massRatio)));
      } else {
        return std::make_pair(std::shared_ptr<OneLevelParticleMapper<T>>(nullptr),
                              std::make_shared<OneLevelParticleMapper<T>>(this->pGrid));
      }
    }

    MapPtrType superOrSubSampleDM(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {
      if (pGrid->pointsToAnyGrid(toGrids)) {
        GridPtrType newGrid;
        if (super)
          newGrid = std::make_shared<SuperSampleGrid<T>>(this->pGrid, ratio);
        else
          newGrid = std::make_shared<SubSampleGrid<T>>(this->pGrid, ratio);
        return std::make_shared<OneLevelParticleMapper<T>>(newGrid);
      } else {
        return std::make_shared<OneLevelParticleMapper<T>>(this->pGrid);
      }
    }


  };

  template<typename T>
  class TwoLevelParticleMapper : public ParticleMapper<T> {
  private:

    using MapType = ParticleMapper<T>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;
    using typename MapType::ConstGridPtrType;

    MapPtrType pLevel1;
    MapPtrType pLevel2;

    GridPtrType pGrid1;
    GridPtrType pGrid2;

    unsigned int n_hr_per_lr;

    size_t totalParticles;
    size_t firstLevel2Particle;

    bool skipLevel1;

    void appendMappedIdsToVector(size_t id0, std::vector<size_t> &ids) const {
      // Finds all the particles at level2 corresponding to the single
      // particle at level1
      T x0, y0, z0;
      std::tie(x0, y0, z0) = pGrid1->getCellCentroid(id0);
      pGrid2->appendIdsInCubeToVector(x0, y0, z0, pGrid1->dx, ids);
    }

    std::vector<size_t> getMappedIds(size_t id0) const {
      std::vector<size_t> ids;
      appendMappedIdsToVector(id0, ids);
      return ids;
    }

    size_t reverseMapId(size_t id2) const {
      auto coord = pGrid2->getCellCentroid(id2);
      return pGrid1->getClosestIdNoWrap(coord);
    }


  public:
    std::vector<size_t> zoomParticleArray; //< the particles on the coarse grid which we wish to replace with their zooms
  protected:
    mutable std::vector<size_t> zoomParticleArrayHiresUnsorted; //< the particles on the fine grid that are therefore included

    void syncL2iterator(const size_t &next_zoom, iterator &level2iterator) const {

      if (zoomParticleArrayHiresUnsorted[next_zoom] > level2iterator.i) {
        level2iterator += zoomParticleArrayHiresUnsorted[next_zoom] - level2iterator.i;
      } else {
        level2iterator -= level2iterator.i - zoomParticleArrayHiresUnsorted[next_zoom];
      }

      assert(level2iterator.i == zoomParticleArrayHiresUnsorted[next_zoom]);
    }

    void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const override {
      indent(s, n);
      s << "i=" << pIterator->i << " into iterator for TwoLevelParticleMapper " << endl;
      iterator &level1iterator = *(pIterator->subIterators[0]);
      iterator &level2iterator = *(pIterator->subIterators[1]);
      const size_t &i = pIterator->i;
      const size_t &next_zoom = pIterator->extraData[0];
      const size_t &next_zoom_index = pIterator->extraData[1];

      if (i >= firstLevel2Particle) {
        indent(s, n);
        s << "Inside level 2 particles" << endl;
        level2iterator.debugInfo(s, n + 1);
      } else {
        indent(s, n);
        s << "Inside level 1 particles" << endl;
        indent(s, n);
        s << "next_zoom = " << next_zoom << "; next_zoom_index = " << next_zoom_index << endl;
        level1iterator.debugInfo(s, n + 1);
      }
    }

    virtual void incrementIterator(iterator *pIterator) const override {

      // set up some helpful shortcut references
      auto &extraData = pIterator->extraData;
      iterator &level1iterator = *(pIterator->subIterators[0]);
      iterator &level2iterator = *(pIterator->subIterators[1]);
      size_t &i = pIterator->i;

      size_t &next_zoom = extraData[0];
      size_t &next_zoom_index = extraData[1];

      // increment the boring-old-counter!
      ++i;

      // now work out what it actually points to...

      if (i >= firstLevel2Particle) {
        // do zoom particles

        if (i == firstLevel2Particle)
          next_zoom = 0;
        else
          ++next_zoom;

        syncL2iterator(next_zoom, level2iterator);


      } else {
        // do normal particles
        //
        // N.B. should never reach here if level1iterator is NULL

        assert((&level1iterator) != nullptr);

        ++level1iterator;

        while (level1iterator.i == next_zoom_index) {
          // we DON'T want to return this particle.... it will be 'zoomed'
          // later on... ignore it
          ++level1iterator;

          // load the next zoom particle
          ++next_zoom;
          if (next_zoom < zoomParticleArray.size())
            next_zoom_index = zoomParticleArray[next_zoom];
          else
            next_zoom_index = size() + 1; // i.e. there isn't a next zoom index!

          // N.B. now it's possible we our 'next' particle is also to be zoomed on,
          // so the while loop goes back round to account for this
        }


      }

    }


    virtual void incrementIteratorBy(iterator *pIterator, size_t increment) const {
      // could be optimized:
      for (size_t i = 0; i < increment; ++i)
        incrementIterator(pIterator);

    }

    virtual void dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
      if (pIterator->i >= firstLevel2Particle)
        pIterator->subIterators[1]->deReference(gp, i);
      else
        pIterator->subIterators[0]->deReference(gp, i);
    }

    void calculateHiresParticleList() const {

      zoomParticleArrayHiresUnsorted.reserve(zoomParticleArray.size() * n_hr_per_lr);

      for (size_t i = 0; i < zoomParticleArray.size(); ++i) {
        appendMappedIdsToVector(zoomParticleArray[i], zoomParticleArrayHiresUnsorted);

      }

      if (!pLevel2->supportsReverseIterator()) {
        // underlying map can't cope with particles being out of order - sort them
        std::sort(zoomParticleArrayHiresUnsorted.begin(), zoomParticleArrayHiresUnsorted.end());
      }

    }


  public:

    virtual iterator begin(const AbstractMultiLevelParticleGenerator<T> &generator) const override {
      if (zoomParticleArray.size() == 0)
        return this->end(generator);

      iterator x(this, generator);
      x.extraData.push_back(0); // current position in zoom ID list
      x.extraData.push_back(zoomParticleArray[0]); // next entry in zoom ID list

      if (!skipLevel1)
        x.subIterators.emplace_back(new iterator(pLevel1->begin(generator)));
      else
        x.subIterators.emplace_back(nullptr);

      x.subIterators.emplace_back(new iterator(pLevel2->begin(generator)));

      syncL2iterator(0, *(x.subIterators.back()));

      return x;
    }


    TwoLevelParticleMapper(MapPtrType &pLevel1,
                           MapPtrType &pLevel2,
                           const std::vector<size_t> &zoomParticles,
                           int n_hr_per_lr,
                           bool skipLevel1 = false) :
      pLevel1(pLevel1), pLevel2(pLevel2),
      pGrid1(pLevel1->getCoarsestGrid()),
      pGrid2(pLevel2->getCoarsestGrid()),
      n_hr_per_lr(n_hr_per_lr),
      skipLevel1(skipLevel1),
      zoomParticleArray(zoomParticles) {

      assert(zoomParticleArray.size() > 0);

      std::sort(zoomParticleArray.begin(), zoomParticleArray.end());

      totalParticles = pLevel1->size() + (n_hr_per_lr - 1) * zoomParticleArray.size();

      firstLevel2Particle = pLevel1->size() - zoomParticleArray.size();

      if (skipLevel1) {
        firstLevel2Particle = 0;
        totalParticles = n_hr_per_lr * zoomParticleArray.size();
      }

      assert(pLevel1->size_gas() == 0);
      assert(pLevel2->size_gas() == 0);

      calculateHiresParticleList();


    }

    virtual void debugInfo(std::ostream &s, int level = 0) const override {
      indent(s, level);
      s << "TwoLevelParticleMapper, n_hr_per_lr=" << n_hr_per_lr << ", firstLevel2Particle=" << firstLevel2Particle <<
        std::endl;
      indent(s, level);
      s << "                      , zoom.size=" << zoomParticleArray.size() << ", zoomed.size=" <<
        zoomParticleArrayHiresUnsorted.size() << std::endl;
      if (skipLevel1) {
        indent(s, level);
        s << "low-res part will be skipped but notionally is:" << endl;
      }
      pLevel1->debugInfo(s, level + 1);
      indent(s, level);
      s << "TwoLevelParticleMapper continues with high-res particles:" << std::endl;
      pLevel2->debugInfo(s, level + 1);
      indent(s, level);
      s << "TwoLevelParticleMapper ends" << std::endl;
    }

    bool references(GridPtrType grid) const override {
      return pLevel1->references(grid) || pLevel2->references(grid);
    }

    virtual void clearParticleList() override {
      pLevel1->clearParticleList();
      pLevel2->clearParticleList();
    }

    void flagParticles(const std::vector<size_t> &genericParticleArray) override {

      std::vector<size_t> grid1particles;
      std::vector<size_t> grid2particles;

      if (pLevel1->size() < zoomParticleArray.size())
        throw std::runtime_error("Zoom particle list is longer than the grid it refers to");

      size_t i0 = pLevel1->size() - zoomParticleArray.size();
      size_t lr_particle;

      size_t last_val = 0;
      size_t numSkippedParticlesLR = 0;
      size_t nZoomParticles = zoomParticleArray.size();

      for (auto i = genericParticleArray.begin(); i != genericParticleArray.end(); ++i) {
        if (*i < last_val)
          throw std::runtime_error("The particle list must be in ascending order");
        last_val = (*i);

        if ((*i) < i0) {
          // particle in low res region. We need to count how many skipped
          // particles there are to work out the address in the original
          // file ordering
          while (numSkippedParticlesLR < nZoomParticles &&
                 zoomParticleArray[numSkippedParticlesLR] <= (*i) + numSkippedParticlesLR)
            ++numSkippedParticlesLR;

          grid1particles.push_back(*i + numSkippedParticlesLR);

        } else {
          // particle in high res region

          if (((*i) - i0) / n_hr_per_lr >= zoomParticleArray.size())
            throw std::runtime_error("Particle ID out of range");

          // find the low-res particle. Note that we push it onto the list
          // even if it's already on there to get the weighting right and
          // prevent the need for a costly search
          lr_particle = zoomParticleArray[((*i) - i0) / n_hr_per_lr];

          grid1particles.push_back(lr_particle);

          // get all the HR particles
          std::vector<size_t> hr_particles = getMappedIds(lr_particle);
          assert(hr_particles.size() == n_hr_per_lr);

          // work out which of these this particle must be and push it onto
          // the HR list
          int offset = ((*i) - i0) % n_hr_per_lr;
          grid2particles.push_back(hr_particles[offset]);
        }

      }

      // Some routines need these lists to be sorted (e.g. getFlaggedParticles below,
      // or in principle there could be an underlying further processing step.
      // So it's time to do it...
      std::sort(grid1particles.begin(), grid1particles.end());
      std::sort(grid2particles.begin(), grid2particles.end());

      // Get the grids to interpret their own particles. This almost
      // certainly just involves making a copy of our list (notwithstanding
      // point made above about a possible underlying processing step). In fact, we're
      // done with our list so they might as well steal the data.
      pLevel1->flagParticles(std::move(grid1particles));
      pLevel2->flagParticles(std::move(grid2particles));

    }


    void getFlaggedParticles(std::vector<size_t> &particleArray) const override {

      // translate level 1 particles - just need to exclude the zoomed particles
      std::vector<size_t> grid1particles;
      pLevel1->getFlaggedParticles(grid1particles);

      size_t zoom_i = 0;
      size_t len_zoom = zoomParticleArray.size();
      for (const size_t &i_lr : grid1particles) {
        while (zoom_i < len_zoom && zoomParticleArray[zoom_i] < i_lr)
          ++zoom_i;

        if (zoom_i == len_zoom || zoomParticleArray[zoom_i] != i_lr) {
          // not a zoom particle: record it in the low res region
          particleArray.push_back(i_lr - zoom_i);
        }
      }

      std::vector<size_t> grid2particles;
      pLevel2->getFlaggedParticles(grid2particles);

      // get the ordering of the level 2 particles
      std::vector<size_t> sortIndex = argsort(zoomParticleArrayHiresUnsorted);


      size_t zoomed_i = 0;
      size_t len_zoomed = zoomParticleArrayHiresUnsorted.size();
      for (const size_t &i_hr : grid2particles) {
        // find the i_hr in the zoomed particle list
        while (zoomed_i < len_zoomed && zoomParticleArrayHiresUnsorted[sortIndex[zoomed_i]] < i_hr)
          ++zoomed_i;

        // If the marked particle is not actually in the output list, ignore it.
        //
        // Older versions of the code throw an exception instead
        if (zoomed_i == len_zoomed || zoomParticleArrayHiresUnsorted[sortIndex[zoomed_i]] != i_hr)
          continue;

        particleArray.push_back(sortIndex[zoomed_i] + firstLevel2Particle);
      }

      std::sort(particleArray.begin(), particleArray.end());

    }

    virtual GridPtrType getCoarsestGrid() override {
      return pLevel1->getCoarsestGrid();
    }

    GridPtrType getFinestGrid() override {
      return pLevel2->getFinestGrid();
    }

    virtual size_t size() const override {
      return totalParticles;
    }

    std::pair<MapPtrType, MapPtrType> addGas(T massRatio, const std::vector<GridPtrType> &toGrids) override {
      bool newskip = skipLevel1;

      auto newLevel1 = pLevel1->addGas(massRatio, toGrids);
      auto newLevel2 = pLevel2->addGas(massRatio, toGrids);

      auto gasSubLevel1 = newLevel1.first;
      auto gasSubLevel2 = newLevel2.first;
      auto dmSubLevel1 = newLevel1.second;
      auto dmSubLevel2 = newLevel2.second;

      decltype(gasSubLevel1) newGasMap;
      decltype(gasSubLevel1) newDmMap;

      if (gasSubLevel1 == nullptr) {
        gasSubLevel1 = pLevel1;
        newskip = true;
      }

      if (gasSubLevel2 != nullptr)
        newGasMap = std::make_shared<TwoLevelParticleMapper<T>>(
          gasSubLevel1, gasSubLevel2, zoomParticleArray,
          n_hr_per_lr, newskip);
      else
        newGasMap = nullptr;

      newDmMap = std::make_shared<TwoLevelParticleMapper<T>>(
        dmSubLevel1, dmSubLevel2, zoomParticleArray,
        n_hr_per_lr, skipLevel1);

      return std::make_pair(newGasMap, newDmMap);


    }

    MapPtrType superOrSubSampleDM(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {

      auto ssub1 = pLevel1->superOrSubSampleDM(ratio, toGrids, super);
      auto ssub2 = pLevel2->superOrSubSampleDM(ratio, toGrids, super);

      // Work out the new list of particles to zoom on
      decltype(zoomParticleArray) newZoomParticles;
      ssub1->clearParticleList();
      pLevel1->flagParticles(zoomParticleArray);
      ssub1->getFlaggedParticles(newZoomParticles);

      size_t new_n_hr_per_lr = n_hr_per_lr;

      // Work out the new number of high res particles per low res
      if (super && ssub2.get() != pLevel2.get()) {
        // super-sampling changed high-res grid
        new_n_hr_per_lr = n_hr_per_lr * ratio * ratio * ratio;
      } else if (!super && ssub1.get() != pLevel1.get()) {
        // sub-sampling changed low-res grid
        new_n_hr_per_lr = n_hr_per_lr * (ratio * ratio * ratio);
      } else if (ssub1.get() != pLevel1.get() || ssub2.get() != pLevel2.get()) {
        // For more general changes, n_hr_per_lr might not be calculable.
        // In the most general map, n_hr_per_lr might vary from particle to particle
        // (e.g. if we start allowing multi-level zooms). This would require
        // re-implementation of the entire mapping class to remove the explicit
        // requirement to know n_hr_per_lr. Work for the future...
        throw runtime_error("Supersampling/subsampling error. The change is in a regime which can't be calculated.");
      }


      return std::make_shared<TwoLevelParticleMapper<T>>(
        ssub1, ssub2, newZoomParticles,
        new_n_hr_per_lr,
        skipLevel1);
    }

  };


  template<typename T>
  class AddGasMapper : public ParticleMapper<T> {
  public:
    using MapType = ParticleMapper<T>;
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


    virtual void dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
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
      indent(s, level);
      s << "AddGasMapper";
      if (gasFirst)
        s << ", gas first:" << std::endl;
      else
        s << ", DM first:" << std::endl;

      firstMap->debugInfo(s, level + 1);

      indent(s, level);

      s << "AddGasMapper";
      if (gasFirst)
        s << ", DM second:" << std::endl;
      else
        s << ", gas second:" << std::endl;

      secondMap->debugInfo(s, level + 1);

      indent(s, level);
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
      firstMap(pFirst), secondMap(pSecond), gasFirst(gasFirst), nFirst(pFirst->size()), nSecond(pSecond->size()) {
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

    virtual void clearParticleList() override {
      firstMap->clearParticleList();
      secondMap->clearParticleList();
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

    virtual iterator begin(const AbstractMultiLevelParticleGenerator<T> &generator) const override {
      iterator i(this, generator);
      i.subIterators.emplace_back(new iterator(firstMap->begin(generator)));
      i.subIterators.emplace_back(new iterator(secondMap->begin(generator)));
      return i;
    }

    virtual iterator beginDm(const AbstractMultiLevelParticleGenerator<T> &generator) const override {
      return (gasFirst ? secondMap : firstMap)->begin(generator);
    }

    virtual iterator endDm(const AbstractMultiLevelParticleGenerator<T> &generator) const override {
      return (gasFirst ? secondMap : firstMap)->end(generator);
    }

    virtual iterator beginGas(const AbstractMultiLevelParticleGenerator<T> &generator) const override {
      return (gasFirst ? firstMap : secondMap)->begin(generator);
    }

    virtual iterator endGas(const AbstractMultiLevelParticleGenerator<T> &generator) const override {
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

      return std::make_shared<AddGasMapper<T>>(
        ssub1, ssub2, gasFirst);
    }


  };

  template
  class MapperIterator<double>;
}

#endif // __MAPPER_HPP
