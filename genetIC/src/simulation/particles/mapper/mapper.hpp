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

				*/
        template<typename GridDataType>
        class ParticleMapper {
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

            virtual void incrementIterator(iterator *pIterator) const {
                ++(pIterator->i);
            }

            virtual void incrementIteratorBy(iterator *pIterator, size_t increment) const {
                pIterator->i += increment;
            }

            virtual void decrementIteratorBy(iterator*, size_t) const {
                throw std::runtime_error("Attempting to reverse in a mapper that does not support random access");
            }

            virtual void dereferenceIterator(const iterator*, ConstGridPtrType &, size_t &) const {
                throw std::runtime_error("There is no grid associated with this particle mapper");
            }

        public:

            virtual void debugInfo(std::ostream &s, int level = 0) const {
                tools::indent(s, level);
                s << "Abstract MapperIterator (shouldn't be here)" << std::endl;
            }

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

            friend std::ostream &operator<<(std::ostream &stream, const ParticleMapper<GridDataType> &I) {
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

            virtual void unflagAllParticles() {

            }

            virtual void flagParticles(const std::vector<size_t>&) {
                throw std::runtime_error("Cannot interpret particles yet; no particle->grid mapper available");
            }

            virtual void flagParticles(std::vector<size_t> &&genericParticleArray) {
                // Sometimes it's helpful to have move semantics, but in general just call
                // the normal implementation
                this->flagParticles(genericParticleArray);
            }

            virtual void getFlaggedParticles(std::vector<size_t>&) const {
                throw std::runtime_error("Cannot get particles; no particle->grid mapper available");
            }

            virtual void extendParticleListToUnreferencedGrids(
                    multilevelcontext::MultiLevelContextInformation<GridDataType> &grids) {
                /* For any grid that is _not_ referenced by this mapper, generate cell flags by matching
                 * to the finest level available in this mapper.
                 *
                 * This is used when specifying constraints wrt an unzoomed simulation in a
                 * zoomed simulation - the cell flags will not reach onto the finest level until this routine
                 * is run.
                 */
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

            virtual iterator begin(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
                return iterator(this, generator);
            }

            virtual bool references(GridPtrType grid) const {
                return getCoarsestGrid().get() == grid.get();
            }

            virtual iterator end(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
                iterator x(this, generator);
                x.i = size();

                return x;
            }

            virtual iterator beginDm(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
                return begin(generator);
            }

            virtual iterator endDm(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
                return end(generator);
            }

            virtual iterator beginGas(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
                return end(generator);
            }

            virtual iterator endGas(const AbstractMultiLevelParticleGenerator<GridDataType> &generator) const {
                return end(generator);
            }


            virtual bool supportsReverseIterator() {
                return false;
            }

            virtual std::pair<MapPtrType, MapPtrType> addGas(T /*massratio*/, const std::vector<GridPtrType>& /*&toGrids*/) {
                throw std::runtime_error("Don't know how to add gas in this context");
            }

            virtual MapPtrType
            superOrSubSampleDM(int /*ratio*/, const std::vector<GridPtrType>& /*&toGrids*/, bool /*super = true*/) {
                throw std::runtime_error("Don't know how to supersample in this context");
            }


        };
    }
};
#endif // __MAPPER_HPP
