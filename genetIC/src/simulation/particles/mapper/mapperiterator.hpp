#ifndef IC_MAPPER_IMPLEMENTATIONS_HPP
#define IC_MAPPER_IMPLEMENTATIONS_HPP

#include <memory>
#include <typeinfo>
#include <src/simulation/particles/generator.hpp>

namespace particle {
  template<typename GridDataType>
  class AbstractMultiLevelParticleGenerator;

  template<typename GT>
  class ParticleGenerator;

  namespace mapper {
    using std::endl;
    using std::cerr;

    // forward declarations:
    template<typename GT>
    class ParticleMapper;

    template<typename GT>
    class OneLevelParticleMapper;

    template<typename GT>
    class TwoLevelParticleMapper;

    template<typename GT>
    class AddGasMapper;

    /*!
     \class MapperIterator
     \brief Class helping the mapper keeping track
     of where within a list of particles a current process (e.g. writing a file) has got to

     Iterators store their current position, but a particle mapper is required to interpret it
     and link the iterator to a cell on the grid. Iterators know about a single grid, but they
     can also contain sub-iterators that refer to other grids in the multi-level context. Each
     iterator is also tied to a specific particle generator, which can be used to quickly
     generate particles by iterating through the structure of grid.

    */
    template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
    class MapperIterator {
    protected:
      // This is a generic structure used by all subclasses to keep
      // track of where they are in a sequential operation

      using MapType = ParticleMapper<GridDataType>;
      using MapPtrType = std::shared_ptr<MapType>;
      using GridType = grids::Grid<T>;
      using ConstGridPtrType = std::shared_ptr<const grids::Grid<T>>;
      using DereferenceType = std::pair<ConstGridPtrType, size_t>;

      friend class ParticleMapper<GridDataType>;

      friend class OneLevelParticleMapper<GridDataType>;

      friend class TwoLevelParticleMapper<GridDataType>;

      friend class AddGasMapper<GridDataType>;


      size_t i; //!< Current position of the iterator.
      std::vector<std::shared_ptr<MapperIterator<GridDataType>>> subIterators; //!< For iterators that have multiple levels, or point to dark matter and baryons separately, stores sub iterators that will be incremented
      std::vector<size_t> extraData; //!< Extra data that may need to be stored by the iterator

      const ParticleMapper<GridDataType> *pMapper; //!< Pointer to the particle mapper that uses the iterator
      const AbstractMultiLevelParticleGenerator<GridDataType> &generator; //!< Generator used by the iterator to create particles for a given cell
      mutable ConstGridPtrType pLastGrid; //!< Pointer to the last grid pointed to
      mutable std::shared_ptr<const ParticleEvaluator<GridDataType>> pLastGridEvaluator; //!< Evaluator for fields on the last grid pointed to

      /*! \brief Constructor, that accepts a particle mapper, and a particle generator
        \param pMapper - particle mapper that uses this iterator
        \param generator - particle generator that will be used to create particles at each cell
      */
      MapperIterator(const ParticleMapper<GridDataType> *pMapper,
                     const AbstractMultiLevelParticleGenerator<GridDataType> &generator) :
        i(0), pMapper(pMapper), generator(generator) {}


    public:

      //! Constructor that copies another iterator
      MapperIterator(const MapperIterator<GridDataType> &source) :
        i(source.i), extraData(source.extraData), pMapper(source.pMapper),
        generator(source.generator), pLastGrid(nullptr) {
        for (const auto &subIterator: source.subIterators) {
          if (subIterator == nullptr)
            subIterators.push_back(nullptr);
          else
            subIterators.push_back(std::make_shared<MapperIterator<GridDataType>>(*subIterator));
        }
      }

      //! Increments the iterator by one step
      MapperIterator &operator++() {
        pMapper->incrementIterator(this);
        return (*this);
      }

      //! Increments the iterator by the specified number of steps
      MapperIterator operator++(int) {
        MapperIterator R = (*this);
        pMapper->incrementIterator(this);
        return R;
      }

      //! Increments the iterator by the specified number of steps
      MapperIterator &operator+=(size_t m) {
        pMapper->incrementIteratorBy(this, m);
        return (*this);
      }

      //! Decrements the iterator by the specified number of steps
      MapperIterator &operator-=(size_t m) {
        pMapper->decrementIteratorBy(this, m);
        return (*this);
      }

      //! Dereferences the iterator at its current position and returns a pointer to the level current pointed at, and the index of the cell pointed to on that level
      DereferenceType operator*() const {
        ConstGridPtrType gp;
        size_t i;
        deReference(gp, i);
        return std::make_pair(gp, i);
      }


      /*! \brief Performs the dereferencing operation
        \param gp - reference to store pointer to level the iterator is pointing to
        \param id - reference to store index for cell pointed to on the relevant level
      */
      void deReference(ConstGridPtrType &gp, size_t &id) const {
        // TODO: weird that this is now partly updating internal state and partly returning results
        pMapper->dereferenceIterator(this, gp, id);
        if (gp != pLastGrid) {
          pLastGrid = gp;
          updateGridReference();
        }
      }


      //! Returns the particle pointed to by the iterator
      Particle<T> getParticle() const {
        ConstGridPtrType pGrid;
        size_t id;
        deReference(pGrid, id);
        return pLastGridEvaluator->getParticle(id);
      }

    protected:
      //! Updates the stored evaluator to match the last grid pointed to.
      void updateGridReference() const {
        pLastGridEvaluator = generator.makeParticleEvaluatorForGrid(*pLastGrid);
      }

      mutable const fields::MultiLevelField<GridDataType> *lastMLField; //!< Pointer to the last multi-level field pointed to
      mutable ConstGridPtrType lastGridPtr; //!< Pointer to the last grid pointed to
      mutable std::shared_ptr<fields::EvaluatorBase<GridDataType, T>> lastEvaluator; //!< Pointer to the last evaluator pointed to.


      /*! \brief Returns the evaluator for the specified multi-level field and grid
        \param multiLevelField - multi-level field to get evaluator for
        \param gridPtr - grid to get evaluator for
      */
      std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>
      getEvaluatorForFieldAndGrid(const fields::MultiLevelField<GridDataType> &multiLevelField,
                                  ConstGridPtrType gridPtr) const {
        // TODO: this evaluator requries optimisation and neatining. In particular it will be extremely inefficient
        // if more than one field is being evaluated at each iteration, because it will need to keep calling
        // makeEvaluator.
        if (lastMLField != &multiLevelField || gridPtr != lastGridPtr) {
          lastEvaluator = fields::makeEvaluator(multiLevelField, *gridPtr);
          lastMLField = &multiLevelField;
          lastGridPtr = gridPtr;
        }
        return lastEvaluator;
      }

    public:

      //! Evaluates the specified multi-level field at the cell currently pointed to by the iterator
      template<typename S>
      auto getField(const fields::MultiLevelField<S> &multiLevelField) const {
        ConstGridPtrType grid_ptr;
        size_t grid_index;
        std::tie(grid_ptr, grid_index) = **this;
        auto evaluator = getEvaluatorForFieldAndGrid(multiLevelField, grid_ptr);

        return (*evaluator)[grid_index];
      }

      //! Returns the index to which the iterator currently points.
      size_t getIndex() const {
        return i;
      }

      //! Returns the number of particles remaining before the end of the particle list, from the current position
      size_t getNumRemainingParticles() const {
        return pMapper->size() - getIndex();
      }

      //! Iterates in parallel, applying the callback function
      size_t parallelIterate(std::function<void(size_t, const MapperIterator &)> callback, size_t nMax) {
        if (pMapper == nullptr) return 0;

        size_t n = std::min(pMapper->size() - i, nMax);
        size_t final_i = i + n;

        if (n == 0) return 0;

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
          if ((unsigned) thread_num < n) {
            (*pThreadLocalIterator) += thread_num;

            for (size_t local_i = thread_num; local_i < n; local_i += num_threads) {
              callback(local_i, *pThreadLocalIterator);
              if (local_i + num_threads < n) (*pThreadLocalIterator) += num_threads;
            }
          }

          if (thread_num != 0)
            delete pThreadLocalIterator;

        }
        if (final_i > i)
          (*this) += final_i - i;
        return n;
      }


      //! Gets the mass in the cell currently pointed at.
      T getMass() const {
        ConstGridPtrType pGrid;
        size_t id;
        deReference(pGrid, id);
        return pLastGridEvaluator->getMass();
      }

      //! Returns a pointer to the pair obtained by dereferencing the iterator at its current position
      std::unique_ptr<DereferenceType> operator->() const {
        return std::unique_ptr<DereferenceType>(new DereferenceType(**this));
      }

      //! Returns true if the two iterators are at the same point, and have the same mapper
      friend bool operator==(const MapperIterator<GridDataType> &lhs, const MapperIterator<GridDataType> &rhs) {
        return (lhs.i == rhs.i) && (lhs.pMapper == rhs.pMapper);
      }

      //! Returns true if the two iterators are not at the same point, or don't have the same mapper
      friend bool operator!=(const MapperIterator<GridDataType> &lhs, const MapperIterator<GridDataType> &rhs) {
        return (lhs.i != rhs.i) || (lhs.pMapper != rhs.pMapper);
      }

      //! Outputs debug information about the iterator to the specified stream
      void debugInfo(std::ostream &s, int n = 0) const {
        pMapper->debugInfoForIterator(s, n, this);
      }


    };

  }
}


#endif
