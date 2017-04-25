#ifndef IC_MAPPER_IMPLEMENTATIONS_HPP
#define IC_MAPPER_IMPLEMENTATIONS_HPP

#include <memory>
#include <typeinfo>

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


      size_t i;
      std::vector<std::shared_ptr<MapperIterator<GridDataType>>> subIterators;
      std::vector<size_t> extraData;
      const ParticleMapper<GridDataType> *pMapper;
      const AbstractMultiLevelParticleGenerator<GridDataType> &generator;
      mutable ConstGridPtrType pLastGrid;
      mutable std::shared_ptr<const ParticleGenerator<GridDataType>> pLastGridGenerator;

      MapperIterator(const ParticleMapper<GridDataType> *pMapper,
                     const AbstractMultiLevelParticleGenerator<GridDataType> &generator) :
        i(0), pMapper(pMapper), generator(generator) {}


    public:

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
        if (gp != pLastGrid) {
          pLastGrid = gp;
          updateGridReference();
        }
      }


      Particle <T> getParticle() const {
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
      auto getField(const fields::MultiLevelField<S> &multiLevelField) const {
        const auto q = **this;
        return multiLevelField.getFieldForGrid(*q.first)[q.second];
      }


      size_t getNextNParticles(std::vector<Particle < T>>

      &particles) {
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
      size_t parallelIterate(std::function<void(size_t, const MapperIterator &

      )> callback,
      size_t nMax
      ) {
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
	  if(thread_num<n) {
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

      friend bool operator==(const MapperIterator<GridDataType> &lhs, const MapperIterator<GridDataType> &rhs) {
        return (lhs.i == rhs.i) && (lhs.pMapper == rhs.pMapper);
      }

      friend bool operator!=(const MapperIterator<GridDataType> &lhs, const MapperIterator<GridDataType> &rhs) {
        return (lhs.i != rhs.i) || (lhs.pMapper != rhs.pMapper);
      }

      void debugInfo() const {
        debugInfo(cerr, 0);
      }

      void debugInfo(std::ostream &s, int n = 0) const {
        pMapper->debugInfoForIterator(s, n, this);
      }


    };

  }
}


#endif
