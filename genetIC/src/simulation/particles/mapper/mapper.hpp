#ifndef __MAPPER_HPP
#define __MAPPER_HPP

#ifdef _OPENMP

#include <omp.h>

#endif

#include <src/simulation/particles/particle.hpp>
#include "src/simulation/grid/grid.hpp"
#include "src/tools/util_functions.hpp"
#include "src/simulation/field/multilevelfield.hpp"
#include "src/simulation/particles/mapper/mapperiterator.hpp"

namespace particle {
  template<typename GridDataType>
  class AbstractMultiLevelParticleGenerator;


  /*!
  \namespace particle::mapper
  \brief Implements a mapper to keep
  track of the correspondance of particles and grid cells.

   The idea of the classes in this file are to provide a way to map particle IDs
   onto grid locations. This becomes quite complicated when one has multiple
   grids and potentially things like gas particles and so on. By making
   specific classes to do very specific bits of the mapping, hopefully we keep
   this complexity managable while also having a lot of flexibility to handle
   different set-ups.
  */
  namespace mapper {
    using std::endl;
    using std::cerr;


    /*!
     \class ParticleMapper
     \brief Top level interface defining a mapper. Implementations include
     one, two levels and gas mappers.

     Particle mappers actually have two main roles:
     1) - incrementing, decrementing, and moving around iterators on the grid structure
     2) - flagging particles on a grid.

    */
    template<typename GridDataType>
    class ParticleMapper : public std::enable_shared_from_this<ParticleMapper<GridDataType>> {
    public:
      using T = tools::datatypes::strip_complex<GridDataType>;
      using MapType = ParticleMapper<GridDataType>;
      using MapPtrType = std::shared_ptr<MapType>;
      using GridType = grids::Grid<T>;
      using GridPtrType = std::shared_ptr<grids::Grid<T>>;
      using ConstGridPtrType = std::shared_ptr<const grids::Grid<T>>;
      using iterator = MapperIterator<GridDataType>;
      using BaseGeneratorType = AbstractMultiLevelParticleGenerator<GridDataType>;

      friend class MapperIterator<GridDataType>;

    protected:

      //! Increment the specified iterator by one step
      virtual void incrementIterator(iterator *pIterator) const {
        ++(pIterator->i);
      }

      /*! \brief Increment the specified iterator by the specified number of steps
        \param pIterator - iterator to increment
        \param increment - amount to increment pIterator by
      */
      virtual void incrementIteratorBy(iterator *pIterator, size_t increment) const {
        pIterator->i += increment;
      }

      //! Decrement the specified iterator by the specified number of steps. Only implemented by derived mapper classes.
      virtual void decrementIteratorBy(iterator *, size_t) const {
        throw std::runtime_error("Attempting to reverse in a mapper that does not support random access");
      }

      //! Dereference the specified iterator, storing the grid pointer of the level pointed to and index of the cell pointed to on that level. Only implemented by derived mapper classes
      virtual void dereferenceIterator(const iterator *, ConstGridPtrType &, size_t &) const {
        throw std::runtime_error("There is no grid associated with this particle mapper");
      }

      //! Get the particle type of the specified iterator. Implemented only by derived mapper classes
      virtual unsigned int
      gadgetParticleTypeFromIterator(const iterator * /*pIterator*/) const {
        throw std::runtime_error("There is no gadget particle type known for this particle mapper");
      }

    public:

      /*! \brief Outputs debug information about the mapper. Should be over-written by derived classes.
        \param s - stream to output debug information to.
        \param level - level about which to output debug information.
      */
      virtual void debugInfo(std::ostream &s, int level = 0) const {
        tools::indent(s, level);
        s << "Abstract MapperIterator (shouldn't be here)" << std::endl;
      }

      /*! \brief Output debug information about the specified iterator to the specified stream
        \param s - stream to output debug information to
        \param n - Number of times to output "| " to the screen for visual output
        \param pIterator - pointer to iterator about which to output debug information.
      */
      virtual void debugInfoForIterator(std::ostream &s, int n, const iterator *pIterator) const {
        tools::indent(s, n);
        s << "i=" << pIterator->i << " into abstract MapperIterator (\?\?)" << endl;
        for (auto q: pIterator->extraData) {
          tools::indent(s, n);
          s << "data=" << q << endl;
        }
        for (auto q: pIterator->subIterators) {
          tools::indent(s, n);
          if (q != nullptr) {
            s << "subiterator: " << endl;
            q->debugInfo(s, n + 1);
          } else {
            s << "null subiterator" << endl;
          }
        }
      }

      /*! \brief Allows output of debug information directly from a stream
      */
      friend std::ostream &operator<<(std::ostream &stream, const ParticleMapper<GridDataType> &I) {
        I.debugInfo(stream);
        return stream;
      }


      //! Returns the total number of particles mapped to be the mapper
      virtual size_t size() const {
        return 0;
      }

      //! Returns the total number of baryonic particles
      virtual size_t size_gas() const {
        return 0;
      }

      //! Returns the total number of dark matter particles
      virtual size_t size_dm() const {
        return size();
      }

      //! Unflags all the flagged particles that use this mapper (implemented only by derived classes)
      virtual void unflagAllParticles() {

      }

      //! Flag the particles with the given IDs in the mapped view. These IDs must be in ascending order.
      virtual void flagParticles(const std::vector<size_t> &) {
        throw std::runtime_error("Cannot interpret particles yet; no particle->grid mapper available");
      }

      //! Overload, flags particles with the supplied IDs
      virtual void flagParticles(std::vector<size_t> &&genericParticleArray) {
        // Sometimes it's helpful to have move semantics, but in general just call
        // the normal implementation
        this->flagParticles(genericParticleArray);
      }

      //! Copies the IDs of flagged particles to the supplied vector (implemented only by derived classes)
      virtual void getFlaggedParticles(std::vector<size_t> &) const {
        throw std::runtime_error("Cannot get particles; no particle->grid mapper available");
      }

      /*! \brief For any grid that is _not_ referenced by this mapper, generate cell flags by matching to the finest level available in this mapper.
         *
         * This is used when specifying constraints wrt an unzoomed simulation in a
         * zoomed simulation - the cell flags will not reach onto the finest level until this routine
         * is run.
         */
      virtual void extendParticleListToUnreferencedGrids(
        multilevelgrid::MultiLevelGrid<GridDataType> &grids) {

        for (size_t i = 0; i < grids.getNumLevels(); i++) {
          auto pGrid = grids.getGridForLevel(i).shared_from_this();

          if (!this->references(pGrid)) {
            vector<size_t> ar;
            GridPtrType proxyGrid = getFinestGrid()->makeProxyGridToMatch(*pGrid);
            proxyGrid->getFlaggedCells(ar);
            pGrid->flagCells(ar);
          }
        }
      }

      //! Returns a pointer to the coarsest grid associated to this mapper
      virtual GridPtrType getCoarsestGrid() {
        throw std::runtime_error("There is no grid associated with this particle mapper");
      }

      //! Returns a constant pointer to the coarsest grid associated to this mapper
      ConstGridPtrType getCoarsestGrid() const {
        return (const_cast<ParticleMapper *>(this)->getCoarsestGrid());
      }


      //! Returns a pointer to the finest grid associated to this mapper
      virtual GridPtrType getFinestGrid() {
        throw std::runtime_error("There is no grid associated with this particle mapper");
      }

      //! Returns a constant pointer to the finest grid associated to this mapper
      ConstGridPtrType getFinestGrid() const {
        return (const_cast<ParticleMapper *>(this)->getFinestGrid());
      }

      //! Returns an iterator with the specified generator set to the beginning of the particle list
      virtual iterator begin(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
        return iterator(this, generator);
      }

      //! Returns true if the specified grid points to the one associated to this mapper
      virtual bool references(GridPtrType grid) const {
        return getCoarsestGrid()->isProxyFor(grid.get());
      }

      //! Returns an iterator with the specified generator set to the end of the particle list
      virtual iterator end(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
        iterator x(this, generator);
        x.i = size();

        return x;
      }

      /*! \brief Get an iterator for the first particle of the specified gadget type.
       *
       * Can be expensive; you're recommended to instead use iterateParticlesOfType */
      virtual iterator beginParticleType(const AbstractMultiLevelParticleGenerator<GridDataType> & /*generator*/,
                                         unsigned int /*particleType*/) const {
        throw std::runtime_error("There is no gadget particle type associated with this particle mapper");

      }

      /*! \brief Get an iterator for the last particle of the specified gadget type.
       *
       * Can be expensive; you're recommended to instead use iterateParticlesOfType */
      virtual iterator endParticleType(const AbstractMultiLevelParticleGenerator<GridDataType> & /*generator*/,
                                       unsigned int /*particleType*/) const {
        throw std::runtime_error("There is no gadget particle type associated with this particle mapper");

      }

      /*! \brief Iterate over all particles of a specified gadget type.

          \param generator - generator used to create particles.
          \param particle_type - gadget type of the particle to iterate over.
          \param callback - function to be applied to each particle.
      */
      void iterateParticlesOfType(
        const particle::AbstractMultiLevelParticleGenerator<GridDataType> &generator,
        unsigned int particle_type,
        const std::function<void(const iterator &)> &callback) const {
        auto begin = beginParticleType(generator, particle_type);
        auto end = endParticleType(generator, particle_type);
        for (auto i = begin; i != end; ++i) {
          callback(i);
        }
      }

      /*! \brief Iterate over all particles in order of the gadget particle types.
       *
          This may not be the same as the genetIC mapper order.
       * **/
      void iterateParticlesOfAllTypes(
        const std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>>> &generator,
        const std::function<void(const iterator &)> &callback, std::vector<size_t> &gadgetTypeToFieldType) const {
        // This function needs to know how to map generators to particle types - information that is held by redirect. This will be
        // different if, say, we switch baryon transfer functions off, as then everything uses the DM generator.
        for (unsigned int particle_type = 0; particle_type < 6; particle_type++) {
          iterateParticlesOfType(*generator[gadgetTypeToFieldType[particle_type]], particle_type, callback);
        }
      }

      //! Returns an iterator with the specified generator pointing to the beginning of the dark matter
      virtual iterator beginDm(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
        return begin(generator);
      }

      //! Returns an iterator with the specified generator pointing to the end of the dark matter
      virtual iterator endDm(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
        return end(generator);
      }

      //! Returns an iterator with the specified generator pointing to the beginning of the baryons
      virtual iterator beginGas(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
        return end(generator);
      }

      //! Returns an iterator with the specified generator pointing to the end of the baryons
      virtual iterator endGas(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
        return end(generator);
      }


      //! Returns true if the mapper supports iterating in reverse
      virtual bool supportsReverseIterator() {
        return false;
      }

      //! Returns a pair of mappers, with the specified fraction of mass in the first and the remainder in the second.
      //! Only the listed grids are targeted; particles based on other grids are only present in the second mapper returned.
      virtual std::pair<MapPtrType, MapPtrType>
      splitMass(T /*massratio*/, const std::vector<GridPtrType> & /*&toGrids*/) {
        throw std::runtime_error("Don't know how to add gas in this context");
      }

      //! Sub/super samples the dark matter pointed to by the mapper.
      virtual MapPtrType
      superOrSubSample(int /*ratio*/, const std::vector<GridPtrType> & /*&toGrids*/, bool /*super = true*/) {
        throw std::runtime_error("Don't know how to supersample in this context");
      }

      virtual MapPtrType
      insertIntermediateResolutionPadding(size_t /* gridCellRatio */, size_t /* padCells */) {
        throw std::runtime_error("Don't know how to insert intermediate resolution padding in this context");
      }

      //! Decouple the flags on all levels, even if they point to the same logical grid.
      //! This is useful
      virtual MapPtrType withIndependentFlags() = 0;

      //! Couple the flags on levels which are virtual copies of each other, which is the default behaviour.
      virtual MapPtrType withCoupledFlags() = 0;


    };
  }
};
#endif // __MAPPER_HPP
