#ifndef IC_TWOLEVELMAPPER_HPP
#define IC_TWOLEVELMAPPER_HPP

#include "src/simulation/particles/mapper/mapper.hpp"
#include <algorithm>

namespace particle {

  namespace mapper {
    using std::endl;
    using std::cerr;

    /*! \class TwoLevelParticleMapper
       *  \brief The critical class for mapping different grids into an assembled zoom simulation
       *
       * This class points to two other mappers, referred to as level 1 and level 2. The level 2 mapper _must_ be
       * a OneLevelParticleMapper, i.e. it can only point to a single grid. The level 1 mapper on the other hand
       * can itself be a TwoLevelParticleMapper, allowing multi-level maps to be built through a left-expanding
       * tree structure.
       *
       * This means that level 1 particles have both a particle ID and a cell ID which may in general be different,
       * whereas level 2 particles always have a particle and cell ID that match.
       */
    template<typename GridDataType>
    class TwoLevelParticleMapper : public ParticleMapper<GridDataType> {

    private:

      using MapType = ParticleMapper<GridDataType>;
      using typename MapType::T;
      using typename MapType::MapPtrType;
      using typename MapType::GridPtrType;
      using typename MapType::GridType;
      using typename MapType::iterator;
      using typename MapType::ConstGridPtrType;

      MapPtrType pLevel1;
      MapPtrType pLevel2;

      GridPtrType pGrid1;
      GridPtrType pGrid2;

      size_t n_hr_per_lr; //!< Number of level 2 particles per replaced level 1 particle

      size_t totalParticles; //!< Total number of particles in the map
      size_t firstLevel2Particle; //!< Index of the first level 2 particle in the map

      bool skipLevel1; //!< Flag to indicate that the map should not include any particles from level 1, only the particles from level 2

      std::vector<size_t> level1ParticlesToReplace; //!< the particles on the level-1 (finest available) grid, which we wish to replace with their zooms
      std::vector<size_t> level1CellsToReplace; //!< the cells on the level-1 (finest available) grid, which we wish to replace with their zooms
      mutable std::vector<size_t> zoomParticleArrayHiresUnsorted; //< the particles/cells on the level-2 grid that are included

      /*! \brief Write n_hr_per_lr level 2 ids into a vector at the given starting index, corresponding to the level 1 cell id.
         *  \param id - level 1 id
            \param start - iterator for the vector to store level 2 ids in
         *
         * Note that level 2 ids are the same whether we're talking about the level 2 mapper or level 2 grid (because
         * level 2 must be a OneLevelMapper).
         *
         * It is the responsibility of the caller to be sure there is enough space reserved in the vector
         */
      void insertLevel2IdsFromLevel1CellId(size_t id, std::vector<size_t>::iterator start) const {
        T x0, y0, z0;
        std::tie(x0, y0, z0) = pGrid1->getCentroidFromIndex(id);
        pGrid2->insertCubeIdsIntoVector(x0, y0, z0, pGrid1->cellSize, start);
      }

    public:

      /*! \brief Returns an iterator pointing at the end of the list of particles with the specified type
        \param generator - generator to be used to create particles
        \param particleType - gadget type of the particles we want to find the end of

        Note that this function has the effect of going to the end of the particles lists for both
        the levels that the two level mapper points to.
      */
      virtual iterator endParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                       unsigned int particleType) const override {
        return iteratorFromSubiterators(pLevel1->endParticleType(generator, particleType),
                                        pLevel2->endParticleType(generator, particleType));
      }

      /*! \brief Returns an iterator pointing at the beginning of the list of particles with the specified type
        \param generator - generator to be used to create particles
        \param particleType - gadget type of the particles we want to find the beginning of

        Note that this function has the effect of going to the begining of the particles lists for both
        the levels that the two level mapper points to.
      */
      virtual iterator beginParticleType(const AbstractMultiLevelParticleGenerator <GridDataType> &generator,
                                         unsigned int particleType) const override {
        return iteratorFromSubiterators(pLevel1->beginParticleType(generator, particleType),
                                        pLevel2->beginParticleType(generator, particleType));
      }

      /*! \brief Returns an iterator pointing to the beginning of the particle list, for both levels pointed to
        \param generator - particle generator to use for this iterator
      */
      virtual iterator begin(const AbstractMultiLevelParticleGenerator <GridDataType> &generator) const override {
        return iteratorFromSubiterators(pLevel1->begin(generator), pLevel2->begin(generator));
      }

    protected:
      /*! \brief Create an iterator for this mapper, given iterators for the two submappers.
       *  \param l1Iterator - iterator for the first level (a one or two level mapper)
       *  \param l2Iterator - iterator for the second level (must be a one level mapper)
       *
       * First, the level 1 iterator is inspected. If it is not already at the end, we are assumed to be in a state
       * where level 1 particles are being iterated over (regardless of the state of the level 2 iterator). The level
       * 1 iterator is adjusted if it's not already on an unzoomed particle, otherwise it's left alone.  The index
       * of the new iterator is set to the corresponding position. The level 2 iterator is reset to the
       * appropriate starting position.
       *
       * If level 1 iterator is at its endpoint, the default assumption is that we are at the *start* of the zoom
       * particles, i.e. the level 2 iterator is NOT inspected. This is because it is an expensive and, in general,
       * not totally well defined operation to go from a general point in the level 2 iterator to a unique
       * point in the overarching iterator, because no particular ordering is assured.
       *
       * The one exception is that if the level 2 iterator points to a null mapper, we are being requested NOT to
       * return any level 2 particles, so we store a nullptr in place of the level 2 iterator to be sure of this.
       *
       */
      iterator iteratorFromSubiterators(const iterator &l1Iterator, const iterator &l2Iterator) const {
        if (level1ParticlesToReplace.size() == 0)
          return this->end(l1Iterator.generator);

        iterator x(this, l1Iterator.generator);
        x.extraData.push_back(0); // current position in zoom ID list
        x.extraData.push_back(level1ParticlesToReplace[0]); // next entry in zoom ID list

        if (!skipLevel1 && l1Iterator.pMapper != nullptr) {
          x.subIterators.emplace_back(new iterator(l1Iterator));
          moveLevel1SubIteratorPastZoomedParticles(x);
          setIteratorIndexFromLevel1SubIterator(x);
        } else {
          x.i = firstLevel2Particle;
          x.subIterators.emplace_back(nullptr);
        }

        if (l2Iterator.pMapper != nullptr) {
          x.subIterators.emplace_back(new iterator(l2Iterator));
          if (l2Iterator != pLevel2->end(l2Iterator.generator)) {
            // TODO: is the following definitely needed?
            adjustLevel2IteratorForSpecifiedZoomParticle(0, *(x.subIterators.back()));
          } else {
            x.i = size();
          }
        } else {
          x.subIterators.emplace_back(nullptr);
        }

        return x;
      }

      //! Uses the iterators on the deepers levels (level 2) to set the index of iterator x on level 1
      void setIteratorIndexFromLevel1SubIterator(iterator &x) const {
        size_t targetLevel1Index = x.subIterators[0]->i;
        if (targetLevel1Index == 0)
          return;
        size_t numZoomParticlesPriorToTarget =
          std::lower_bound(this->level1ParticlesToReplace.begin(), this->level1ParticlesToReplace.end(),
                           targetLevel1Index)
          - this->level1ParticlesToReplace.begin();

        x.i = targetLevel1Index - numZoomParticlesPriorToTarget;

        if (x.i < firstLevel2Particle) {
          x.extraData[0] = numZoomParticlesPriorToTarget;
        } else {
          x.extraData[0] = 0;
        }
        x.extraData[1] = this->level1ParticlesToReplace[x.extraData[0]];
      }

    public:


      /*! \brief Constructs a two level mapper from two mappers and a vector of particles on the coarse level to replace
        \param pLevel1 - the coarser particle mapper, which can itself have multiple levels
        \param pLevel2 - the fine particle mapper, which can point to only one grid (i.e. must be OneLevelParticleMapper)
        \param level1CellsToReplace - the list of cells on level 1 that will be zoomed. NB these are specified
                                        relative to the finest grid on level 1 (not relative to the level 1 mapper itself)
        \param skipLevel1 - If true, only the level 2 particles are listed, otherwise, level1 particles are listed first, then level 2
      */
      TwoLevelParticleMapper(MapPtrType pLevel1,
                             MapPtrType pLevel2,
                             const std::vector<size_t> &level1CellsToReplace,
                             bool skipLevel1 = false) :
        pLevel1(pLevel1), pLevel2(pLevel2),
        pGrid1(pLevel1->getFinestGrid()),
        pGrid2(pLevel2->getCoarsestGrid()),
        skipLevel1(skipLevel1),
        level1CellsToReplace(level1CellsToReplace) {

        assert(level1CellsToReplace.size() > 0);
        // all level 2 particles must be of the same resolution:
        assert(typeid(*pLevel2) == typeid(OneLevelParticleMapper<GridDataType>));
        assert(pLevel1->size_gas() == 0);
        assert(pLevel2->size_gas() == 0);

        auto pLevel1Indept = pLevel1->withIndependentFlags();

        std::sort(this->level1CellsToReplace.begin(), this->level1CellsToReplace.end());
        pLevel1Indept->unflagAllParticles();
        pLevel1Indept->getFinestGrid()->flagCells(this->level1CellsToReplace);
        pLevel1Indept->getFlaggedParticles(level1ParticlesToReplace);

        /* If the following assertion fails, something is inconsistent. There should be a 1-1 map between the
         * level1CellsToReplace and level1ParticlesToReplace. Possible sources of inconsistency are:
         *
         *  - some of the cells that are flagged are not reachable from the level 1 mapper. For example, if the level 1
         *    mapper is already representing a zoom, only cells in the zoom region should be flagged. This is because
         *    the TwoLevelParticleMapper logic assumes a fixed number of level 2 cells per "missing" level 1 cell.
         *
         *  - the input array level1CellsToReplace is somehow inappropriate.  For example, there was for a long
         *    time a bug where actually the IDs being supplied from various points in the code where level 1 particle
         *    IDs, not cell IDs. This bug was masked in many situations where the level 1 was a OneLevelParticleMapper,
         *    for which there is no distinction between cell and particle IDs.
         * */
        if (this->level1CellsToReplace.size() != level1ParticlesToReplace.size()) {
          logging::entry() << "WANTED:" << endl;
          for (size_t i = 0; i < this->level1CellsToReplace.size(); ++i) {
            logging::entry() << this->level1CellsToReplace[i] << endl;
          }
          pLevel1->unflagAllParticles();
          pLevel1->flagParticles(level1ParticlesToReplace);
          this->level1CellsToReplace.clear();
          pGrid1->getFlaggedCells(this->level1CellsToReplace);
          logging::entry() << "GOT:" << endl;
          for (size_t i = 0; i < this->level1CellsToReplace.size(); ++i) {
            logging::entry() << this->level1CellsToReplace[i] << " (" << this->level1ParticlesToReplace[i] << ")" << endl;
          }
          assert(false);
        }
        n_hr_per_lr = tools::getRatioAndAssertPositiveInteger(pGrid1->cellSize, pGrid2->cellSize);
        n_hr_per_lr *= n_hr_per_lr * n_hr_per_lr;

        totalParticles = pLevel1->size() + (n_hr_per_lr - 1) * level1ParticlesToReplace.size();

        firstLevel2Particle = pLevel1->size() - level1ParticlesToReplace.size();

        if (skipLevel1) {
          firstLevel2Particle = 0;
          totalParticles = n_hr_per_lr * level1ParticlesToReplace.size();
        }

        calculateHiresParticleList();

      }

      /*! \brief Output debug information about this mapper to the specified stream, regarding the specified level
          \param s - stream to output debug information to
          \param level - level about which to output debug information
      */
      virtual void debugInfo(std::ostream &s, int level = 0) const override {
        tools::indent(s, level);
        s << "TwoLevelParticleMapper, n_hr_per_lr=" << n_hr_per_lr << ", firstLevel2Particle="
          << firstLevel2Particle <<
          std::endl;
        tools::indent(s, level);
        s << "                      , zoom.size=" << level1ParticlesToReplace.size() << ", zoomed.size=" <<
          zoomParticleArrayHiresUnsorted.size() << std::endl;
        if (skipLevel1) {
          tools::indent(s, level);
          s << "low-res part will be skipped but notionally is:" << endl;
        }
        pLevel1->debugInfo(s, level + 1);
        tools::indent(s, level);
        s << "TwoLevelParticleMapper continues with high-res particles:" << std::endl;
        pLevel2->debugInfo(s, level + 1);
        tools::indent(s, level);
        s << "TwoLevelParticleMapper ends" << std::endl;
      }

      //! Returns true if either level 1 or level 2 references the specified grid.
      bool references(GridPtrType grid) const override {
        return pLevel1->references(grid) || pLevel2->references(grid);
      }

      //! Unflags all the flagged particles on level 1 and 2
      virtual void unflagAllParticles() override {
        pLevel1->unflagAllParticles();
        pLevel2->unflagAllParticles();
      }

      //! Flags the particles specified by orderedParticleIndices. Flags the level 1 particles first, and then the level 2 particles that lie in the high-res region
      void flagParticles(const std::vector<size_t> &orderedParticleIndices) override {

        std::vector<size_t> level1particles;

        if (pLevel1->size() < level1ParticlesToReplace.size())
          throw std::runtime_error("Zoom particle list is longer than the grid it refers to");

        size_t firstHiresParticleInMapper = pLevel1->size() - level1ParticlesToReplace.size();


        size_t last_val = 0;
        size_t numSkippedParticlesLR = 0;
        size_t nZoomParticles = level1ParticlesToReplace.size();
        size_t firstHrParticleInInput = std::numeric_limits<size_t>::max();


        for (size_t i = 0; i < orderedParticleIndices.size(); ++i) {
          size_t thisParticle = orderedParticleIndices[i];
          if (thisParticle < last_val)
            throw std::runtime_error("The particle list must be in ascending order");
          last_val = (thisParticle);

          if (thisParticle < firstHiresParticleInMapper) {
            // particle in low res region. We need to count how many skipped
            // particles there are to work out the address in the original
            // file ordering. NB this approach assumes the input particle array is
            // in ascending order (and the zoomParticleArray).
            while (numSkippedParticlesLR < nZoomParticles &&
                   level1ParticlesToReplace[numSkippedParticlesLR] <= thisParticle + numSkippedParticlesLR)
              ++numSkippedParticlesLR;

            level1particles.push_back(thisParticle + numSkippedParticlesLR);

          } else {
            // particle in high res region. This is now handled by the separate loop below.
            // The loop for high-res particles is parallelised (unlike the low-res particles)
            // because generating the right indices was found to be slow during profiling.
            firstHrParticleInInput = i;
            break;
          }
        }

        if (firstHrParticleInInput == std::numeric_limits<size_t>::max()) {
          firstHrParticleInInput = orderedParticleIndices.size(); // No HR particles in input
        }


        std::vector<size_t> level2particles;
        level2particles.resize(orderedParticleIndices.size() - firstHrParticleInInput);

        // N.B. the following parallelism does not seem to achieve much speed-up.
        // Is this because it is memory access bound, or some more subtle reason?

#ifdef OPENMP
#pragma omp parallel
        {
#endif
          std::vector<size_t> hrCellsCache;
          std::vector<size_t> localLrParticles;
          hrCellsCache.resize(n_hr_per_lr);
          size_t lrParticleLastAccessed = std::numeric_limits<size_t>::max();

#ifdef OPENMP
#pragma omp for schedule(static, n_hr_per_lr*10)
#endif
          for (size_t i = firstHrParticleInInput; i < orderedParticleIndices.size(); ++i) {
            size_t thisParticle = orderedParticleIndices[i];
            size_t lr_index = (thisParticle - firstHiresParticleInMapper) / n_hr_per_lr;

            if (lr_index >= level1ParticlesToReplace.size())
              throw std::runtime_error("Particle ID out of range");

            // find the low-res particle.

            if (lrParticleLastAccessed != level1ParticlesToReplace[lr_index]) {
              lrParticleLastAccessed = level1ParticlesToReplace[lr_index];
              localLrParticles.push_back(lrParticleLastAccessed);
              // get all the HR particles
              insertLevel2IdsFromLevel1CellId(level1CellsToReplace[lr_index], hrCellsCache.begin());
            }

            // work out which of these this particle must be and push it onto
            // the HR list
            size_t offset = (thisParticle - firstHiresParticleInMapper) % n_hr_per_lr;

            // NB here we assume that level 2 is a OneLevelParticleMapper so that the cell ID can be
            // taken also to be a particle ID. This assumption is explicitly tested with an assert
            // in the constructor.
            level2particles[i - firstHrParticleInInput] = hrCellsCache[offset];

          }

#ifdef OPENMP
#pragma omp critical
#endif
          {
            level1particles.insert(level1particles.end(), localLrParticles.begin(), localLrParticles.end());
          }

#ifdef OPENMP
        }
#endif

        // The divided particle lists will now be passed to the underlying particle mappers,
        // which require the input to be sorted
        std::sort(level1particles.begin(), level1particles.end());
        std::sort(level2particles.begin(), level2particles.end());

        // We're done with our list so they might as well steal the data.
        pLevel1->flagParticles(std::move(level1particles));
        pLevel2->flagParticles(std::move(level2particles));

      }


      //! Copies the ids of flagged particles into particleArray. Guarantees the result is sorted in ascending order.
      void getFlaggedParticles(std::vector<size_t> &particleArray) const override {

        // translate level 1 particles - just need to exclude the zoomed particles
        std::vector<size_t> grid1particles;
        pLevel1->getFlaggedParticles(grid1particles);

        size_t zoom_i = 0;
        size_t len_zoom = level1ParticlesToReplace.size();
        for (const size_t &i_lr : grid1particles) {
          while (zoom_i < len_zoom && level1ParticlesToReplace[zoom_i] < i_lr)
            ++zoom_i;

          if (zoom_i == len_zoom || level1ParticlesToReplace[zoom_i] != i_lr) {
            // not a zoom particle: record it in the low res region
            particleArray.push_back(i_lr - zoom_i);
          }
        }

        std::vector<size_t> grid2particles;
        pLevel2->getFlaggedParticles(grid2particles);

        // get the ordering of the level 2 particles
        std::vector<size_t> sortIndex = tools::argsort(zoomParticleArrayHiresUnsorted);


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

      //! Returns a pointer to the coarsest grid, which is one of the levels of the level 1 grid
      virtual GridPtrType getCoarsestGrid() override {
        return pLevel1->getCoarsestGrid();
      }

      //! Return a pointer to the finest grid, which in this case is the level 2 grid
      GridPtrType getFinestGrid() override {
        // Note - the total mapper always starts with the deepest grid in the stack, with level 1 being the grids above that.
        // This means that getting the finest grid this way from the total mapper will always return the actual finest grid,
        // but if for some reason we called this function on a two-level mapper higher in the stack, it would point only to the
        // deepest grid in that stack, not the deepest in this one. This is different to the getCoarsestGrid function, which
        // always recurses until it reaches the coarsest grid.
        return pLevel2->getFinestGrid();
      }

      //! Returns the total number of particles in the stack
      virtual size_t size() const override {
        return totalParticles;
      }

      /*! \brief Makes two copies of the mapper with differently scaled mass fractions, one for dark matter, one for baryons
          \param massRatio - fraction of total matter mass in the first mapper
          \param toGrids - only adds gas if the grids of the levels pointed to by this mapper correspond to the grids specified here
      */
      std::pair<MapPtrType, MapPtrType> splitMass(T massRatio, const std::vector<GridPtrType> &toGrids) override {
        bool newskip = skipLevel1;

        auto newLevel1 = pLevel1->splitMass(massRatio, toGrids);
        auto newLevel2 = pLevel2->splitMass(massRatio, toGrids);

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
          newGasMap = std::make_shared<TwoLevelParticleMapper<GridDataType>>(
            gasSubLevel1, gasSubLevel2, level1CellsToReplace,
            newskip);
        else
          newGasMap = nullptr;

        newDmMap = std::make_shared<TwoLevelParticleMapper<GridDataType>>(
          dmSubLevel1, dmSubLevel2, level1CellsToReplace,
          skipLevel1);

        return std::make_pair(newGasMap, newDmMap);


      }

      /*! \brief Sub- or super-samples the dark matter
          \param ratio - factor to sub- or super-sample by
          \param toGrids - only perform the sub/super-sampling if the grids match this list
          \param super - super-sample if true, sub-sample if false
      */
      MapPtrType superOrSubSample(int ratio, const std::vector<GridPtrType> &toGrids, bool super) override {

        auto subMapperLR = pLevel1->superOrSubSample(ratio, toGrids, super);
        auto subMapperHR = pLevel2->superOrSubSample(ratio, toGrids, super);

        // Work out the new list of particles to zoom on
        decltype(level1ParticlesToReplace) newLevel1CellsToReplace;
        subMapperLR->unflagAllParticles();
        pLevel1->flagParticles(level1ParticlesToReplace);
        subMapperLR->getFinestGrid()->getFlaggedCells(newLevel1CellsToReplace);

        try {
          auto result = std::make_shared<TwoLevelParticleMapper<GridDataType>>(
            subMapperLR, subMapperHR, newLevel1CellsToReplace,
            skipLevel1);

          return result;
        }
        catch (std::out_of_range &e) {
          // Logically the above operation could fail.
          // See test_15_edge_subsampling for an example of when this might happen
          throw std::runtime_error("Following sub/super-sampling, the zoom particles no longer fit fully on their "
                                   "grids. Try rerunning with strays_on, or (better) increase the size of the zoom "
                                   "region.");
        };
      }

      MapPtrType withIndependentFlags() override {
        return std::make_shared<TwoLevelParticleMapper<T>>(pLevel1->withIndependentFlags(),
                                                           pLevel2->withIndependentFlags(),
                                                           level1CellsToReplace, skipLevel1);
      }

      MapPtrType withCoupledFlags() override {
        return std::make_shared<TwoLevelParticleMapper<T>>(pLevel1->withCoupledFlags(), pLevel2->withCoupledFlags(),
                                                           level1CellsToReplace, skipLevel1);
      }

      MapPtrType insertIntermediateResolutionPadding(size_t ratio, size_t padCells) override {

        auto mapperLR = pLevel1->insertIntermediateResolutionPadding(ratio, padCells);
        auto mapperHR = pLevel2->insertIntermediateResolutionPadding(ratio, padCells);

        assert(mapperHR == pLevel2); // we're expecting level 2 to be a fixed high resolution grid => no changes made

        // level 1 may have been supersampled. We need an updated list of zoom particles relative to the new
        // finest level.
        std::vector<size_t> newLevel1CellsToReplace;
        pLevel1->unflagAllParticles();
        pLevel1->flagParticles(level1ParticlesToReplace);
        mapperLR->getFinestGrid()->getFlaggedCells(newLevel1CellsToReplace);

        // Now, what is the new linear ratio between finest level 1 and level 2
        size_t newCellLengthRatio = tools::getRatioAndAssertInteger(mapperLR->getFinestGrid()->cellSize,
                                                                    mapperHR->getCoarsestGrid()->cellSize);

        if (newCellLengthRatio > ratio) {
          // there is too big a resolution ratio between our levels. Do something about it.
          // The strategy is to supersample the coarsest LR grid, and have a region around our zoom region
          // consisting of such supersampled particles.


          GridPtrType finestLRGrid = mapperLR->getFinestGrid();
          GridPtrType superSampleVersionOfFinestLRGrid = finestLRGrid->makeSupersampled(ratio);

          std::vector<size_t> expandedLevel1CellsToReplace;
          std::vector<size_t> supersampledLevel1CellsToReplace;

          // the cells were already flagged above, when we got newLevel1CellsToReplace
          superSampleVersionOfFinestLRGrid->getFlaggedCells(supersampledLevel1CellsToReplace);

          // Next, we decouple flagging throughout the mapper. This means that flagging on supersampled grids
          // doesn't just flag the whole supercell, which is important to figure out the precise nested virtual
          // zoom levels.
          mapperLR = mapperLR->withIndependentFlags();
          finestLRGrid = mapperLR->getFinestGrid();
          // actually we now want to expand the region -- so that there is space for a 'buffer' region
          finestLRGrid->expandFlaggedRegion(padCells);

          // The simplest thing to do next would be:
          //
          //   finestLRGrid->getFlaggedCells(expandedLevel1CellsToReplace);
          //
          // However this ignores the fact that some cells in the expanded region might not be accessible, i.e.
          // they may be de-refined themselves. Thus, we need to remove any flags on cells that are inaccessible:

          std::vector<size_t> tempList;
          mapperLR->getFlaggedParticles(tempList);
          mapperLR->unflagAllParticles();
          mapperLR->flagParticles(tempList);

          // Finally, we're ready to get our newly expanded region that has the buffer around the original region
          finestLRGrid->getFlaggedCells(expandedLevel1CellsToReplace);
          auto expandedsize = finestLRGrid->getFlaggedCellsPhysicalSize();

          // check we now have a superset of the original flags
          assert(std::includes(expandedLevel1CellsToReplace.begin(), expandedLevel1CellsToReplace.end(),
                               newLevel1CellsToReplace.begin(), newLevel1CellsToReplace.end()));

          // Now we give our new grid independent flags, which is necessary for the recursion below where we may be
          // trying to build a further padded layer around the one we already made. These independent flags will be
          // reverted to coupled flags as a final step below.

          superSampleVersionOfFinestLRGrid = superSampleVersionOfFinestLRGrid->withIndependentFlags();

          /*
           * The following, to make use of the high resolution information where we have it, would be fabulous:
          superSampleVersionOfFinestLRGrid = std::make_shared<grids::ResolutionMatchingGrid<T>>(
            pLevel2->getCoarsestGrid()->makeSubsampled(newCellLengthRatio/ratio),
            finestLRGrid)->withIndependentFlags();

           * However BEWARE. If the simulation has gas, the level 2 grid might be a 'mass scaled grid' whereas the
           * level 1 grid would not be. This can lead to loss of mass conservation, which is very bad indeed.
           *
           * For now, this possibility is disabled.
          */

          // Now ready to assemble our new mapper:
          auto mapperForPaddingRegion = std::make_shared<OneLevelParticleMapper<GridDataType>>(
            superSampleVersionOfFinestLRGrid);
          auto mapperForLowResAndPaddingRegion = std::make_shared<TwoLevelParticleMapper<GridDataType>>(mapperLR,
                                                                                                        mapperForPaddingRegion,
                                                                                                        expandedLevel1CellsToReplace,
                                                                                                        skipLevel1);

          auto mapperWithEverything = std::make_shared<TwoLevelParticleMapper<GridDataType>>(
            mapperForLowResAndPaddingRegion, mapperHR,
            supersampledLevel1CellsToReplace,
            skipLevel1);

          logging::entry() << "Added a padding region of effective resolution " <<
               superSampleVersionOfFinestLRGrid->getEffectiveSimulationSize() << " of physical size ~" <<
               expandedsize << "Mpc/h" << " (inner resolution "
               << mapperHR->getFinestGrid()->getEffectiveSimulationSize() << ", size ~" <<
               superSampleVersionOfFinestLRGrid->getFlaggedCellsPhysicalSize() << "Mpc/h)" << endl;

          // the above might be ready -- but in principle it's possible we need to insert ANOTHER level of padding,
          // which is handled by recursion.
          //
          // Furthermore, we now recouple all the flags as all other parts of the code will expect
          return mapperWithEverything->insertIntermediateResolutionPadding(ratio, padCells)->withCoupledFlags();

        } else {
          // our resolution contrast is fine - generate a new TwoLevelParticleMapper pointing to the new mappers,
          // in case they changed
          return std::make_shared<TwoLevelParticleMapper<GridDataType>>(mapperLR, mapperHR,
                                                                        newLevel1CellsToReplace, skipLevel1);
        }


      }

    protected:


      //!  Outputs debug informaiton about the specified iterator to the specified stream
      void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const override {
        tools::indent(s, n);
        s << "i=" << pIterator->i << " into iterator for TwoLevelParticleMapper " << endl;
        iterator &level1iterator = *(pIterator->subIterators[0]);
        iterator &level2iterator = *(pIterator->subIterators[1]);
        const size_t &i = pIterator->i;
        const size_t &next_zoom = pIterator->extraData[0];
        const size_t &next_zoom_index = pIterator->extraData[1];

        if (i >= firstLevel2Particle) {
          tools::indent(s, n);
          s << "Inside level 2 particles" << endl;
          level2iterator.debugInfo(s, n + 1);
        } else {
          tools::indent(s, n);
          s << "Inside level 1 particles" << endl;
          tools::indent(s, n);
          s << "next_zoom = " << next_zoom << "; next_zoom_index = " << next_zoom_index << endl;
          level1iterator.debugInfo(s, n + 1);
        }
      }

      //! Increments the iterator by incrementing the relevent sub-iterator for the level pointed to
      virtual void incrementIterator(iterator *pIterator) const override {

        // set up some helpful shortcut references
        auto &extraData = pIterator->extraData;
        auto &pLevel1iterator = pIterator->subIterators[0];
        auto &pLevel2iterator = pIterator->subIterators[1];
        size_t &i = pIterator->i;

        size_t &next_zoom = extraData[0];

        // increment the boring-old-counter!
        ++i;

        // now work out what it actually points to...

        if (i >= firstLevel2Particle) {
          // do zoom particles

          if (i == firstLevel2Particle) {
            next_zoom = 0;
          } else {
            ++next_zoom;
          }

          if ((pLevel2iterator) != nullptr)
            adjustLevel2IteratorForSpecifiedZoomParticle(next_zoom, *pLevel2iterator);


        } else {

          ++(*pLevel1iterator);

          moveLevel1SubIteratorPastZoomedParticles(*pIterator);


        }

      }



      // Low-level iterator synchronization methods:

      /*! \brief Moves the specified iterator to the position of the specified particle index (unless it is beyond the end of the list, in which case do nothing)
          \param next_zoom - particle id to move to, defined relative to level 2
          \param level2iterator - iterator for level 2 that needs to be moved
      */
      void adjustLevel2IteratorForSpecifiedZoomParticle(const size_t &next_zoom, iterator &level2iterator) const {

        if (next_zoom >= zoomParticleArrayHiresUnsorted.size()) {
          // beyond end. This is OK because we need to be able to go one beyond the end in an iterator loop
          assert(next_zoom == zoomParticleArrayHiresUnsorted.size());
          return;
        }

        if (zoomParticleArrayHiresUnsorted[next_zoom] > level2iterator.i) {
          level2iterator += zoomParticleArrayHiresUnsorted[next_zoom] - level2iterator.i;
        } else {
          level2iterator -= level2iterator.i - zoomParticleArrayHiresUnsorted[next_zoom];
        }

        assert(level2iterator.i == zoomParticleArrayHiresUnsorted[next_zoom]);
      }

      //! Moves the specified level 1 iterator past the zoom particles on level 2, if we are incremented into that list.
      void moveLevel1SubIteratorPastZoomedParticles(iterator &forIterator) const {
        iterator &level1iterator = *(forIterator.subIterators[0]);
        size_t &next_zoom = forIterator.extraData[0];
        size_t &next_zoom_index = forIterator.extraData[1];
        while (level1iterator.i == next_zoom_index) {
          // we DON'T want to return this particle.... it will be 'zoomed'
          // later on... ignore it
          ++level1iterator;

          // load the next zoom particle
          ++next_zoom;
          if (next_zoom < level1ParticlesToReplace.size())
            next_zoom_index = level1ParticlesToReplace[next_zoom];
          else
            next_zoom_index = size() + 1; // i.e. there isn't a next zoom index!

          // N.B. now it's possible we our 'next' particle is also to be zoomed on,
          // so the while loop goes back round to account for this
        }
      }


      /*! \brief Increments the specified iterator by the specified amount
        \param pIterator - iterator to move
        \param increment - number of steps to increment the iterator by
      */
      virtual void incrementIteratorBy(iterator *pIterator, size_t increment) const override {
        // could be optimized:
        for (size_t i = 0; i < increment; ++i)
          incrementIterator(pIterator);

      }

      /*! \brief Dereference the specified iterator, storing the pointed to grid and cell index
        \param pIterator - iterator to dereference
        \param gp - reference to where the grid pointer should be stored
        \param i - reference to where the cell index should be stored
      */
      virtual void
      dereferenceIterator(const iterator *pIterator, ConstGridPtrType &gp, size_t &i) const override {
        if (pIterator->i >= firstLevel2Particle)
          pIterator->subIterators[1]->deReference(gp, i);
        else
          pIterator->subIterators[0]->deReference(gp, i);
      }


      //! Creates a list of zoom particles from the list of level 1 particles that need to be replaced
      void calculateHiresParticleList() const {
        zoomParticleArrayHiresUnsorted.resize(level1ParticlesToReplace.size() * n_hr_per_lr);

        bool failed = false;

#pragma omp parallel for
        for (size_t i = 0; i < level1CellsToReplace.size(); ++i) {
          try {
            insertLevel2IdsFromLevel1CellId(level1CellsToReplace[i],
                                            zoomParticleArrayHiresUnsorted.begin() + (i * n_hr_per_lr));
          } catch (std::out_of_range &e) {
            failed = true; // OpenMP does not allow exceptions to propagate out of the parallel region :-(
          }
        }

        if (failed) {
          throw std::out_of_range("Requested zoom region falls outside high resolution grid");
        }

        if (!pLevel2->supportsReverseIterator()) {
          // underlying map can't cope with particles being out of order - sort them
          std::sort(zoomParticleArrayHiresUnsorted.begin(), zoomParticleArrayHiresUnsorted.end());
        }

      }


    };
  }
}
#endif
