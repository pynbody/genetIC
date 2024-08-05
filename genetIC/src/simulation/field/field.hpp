#ifndef IC_FIELD_HPP
#define IC_FIELD_HPP

#include <memory>
#include <vector>
#include <cassert>
#include <src/simulation/filters/filter.hpp>
#include <src/tools/numerics/fourier.hpp>
#include <execinfo.h>
#include "src/io/numpy.hpp"
#include "src/simulation/grid/grid.hpp"
#include "src/tools/numerics/tricubic.hpp"
#include "src/tools/lru_cache.hpp"

/*!
    \namespace fields
    \brief Store and manipulate fields of various types, both on individual grids and on multiple grid levels
 */



namespace tools {
  namespace numerics {
    namespace fourier {
      /*! \class FieldFourierManager
          \brief Class for handling Fourier transforms. Has complex and real specialisations.
      */
      template<typename T, typename S=datatypes::strip_complex<T>>
      class FieldFourierManager;
      // implementation in fourier.hpp
      
    }
  }
}


namespace fields {
  namespace cache {

    template<typename T>
    thread_local
      std::unique_ptr<
        tools::lru_cache<std::tuple<int,int,int,const void*>,
                                          numerics::LocalUnitTricubicApproximation<T>
                                          >
                      > cachedInterpolators;
    // The above is a workaround for the fact that thread_local variables are not guaranteed to be initialized
    // (or at least I can't find any documentation that I understand on this).

    thread_local size_t cacheHits, cacheMisses;
    bool enabled;
    size_t accumHits, accumMisses;

    template<typename T>
    void enableInterpolationCaches() {
      enabled = true;
#pragma omp parallel default(none)
      {
        cachedInterpolators<T> = std::make_unique<
          tools::lru_cache<std::tuple<int,int,int,const void*>,
                                      numerics::LocalUnitTricubicApproximation<T>
                                      >
          >(1024);
        cacheHits = 0;
        cacheMisses = 0;
      }
    }

    template <typename T>
    void disableInterpolationCaches() {
      enabled = false;
      cachedInterpolators<T> = nullptr;
#pragma omp parallel default(none) shared(accumHits, accumMisses)
      {
#pragma omp critical
        {
          // pool all hits/misses across threads
          accumHits += cacheHits;
          accumMisses += cacheMisses;
        }

#pragma omp master
        {
#ifdef DEBUG_INFO
          if(accumHits>0 || accumMisses>0) {
            double fracHits = 100 * double(accumHits) / double(accumHits + accumMisses);
            double fracMisses = 100 * double(accumMisses) / double(accumHits + accumMisses);
            logging::entry() << std::setprecision(2) << "Interpolation cache performance report. Hits: " << accumHits
                      << " (" << fracHits << "%); misses: " << accumMisses << " (" << fracMisses << "%)"
                      << std::defaultfloat << std::endl;
          }
#endif
        }
      }
    }
  }

  template<typename D, typename C>
  class EvaluatorBase;

  template<typename D, typename C>
  class Field;

  template<typename D>
  class MultiLevelField;

  // Implementation in evaluator.hpp:
  
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  std::shared_ptr<EvaluatorBase<DataType, CoordinateType>> makeEvaluator(const Field<DataType, CoordinateType> &field,
                                                                         const grids::Grid<CoordinateType> &grid);

  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  std::shared_ptr<EvaluatorBase<DataType, CoordinateType>> makeEvaluator(const MultiLevelField<DataType> &field,
                                                                         const grids::Grid<CoordinateType> &grid);

  size_t peakMemUsage = 0;
  size_t currentMemUsage = 0;

  std::string formatBytes(size_t bytes) {
    std::string suffixes[] = {"B", "KB", "MB", "GB", "TB"};
    int suffix = 0;
    while (bytes > 10*1024 && suffix<4) {
      bytes /= 1024;
      suffix++;
    }
    return std::to_string(bytes) + suffixes[suffix];
  }

  void memUsagePeriodicReportInThread() {
    // Check memory usage every 0.1 second. If zero, report the peak memory usage and return.
    // If non-zero, report the current memory usage but only if that has changed since the
    // last report, and the last report was at least 10 seconds ago

    size_t lastReportedMemUsage = 0;
    auto lastReportTime = std::chrono::steady_clock::now();

    while(true) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      if(currentMemUsage==0) {
        logging::entry() << "Peak memory usage: " << formatBytes(peakMemUsage) << std::endl;
        return;
      }

      if(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - lastReportTime).count() > 10) {
        if(currentMemUsage != lastReportedMemUsage) {
          logging::entry() << "Current memory usage: " << formatBytes(currentMemUsage) << std::endl;
          lastReportedMemUsage = currentMemUsage;
          lastReportTime = std::chrono::steady_clock::now();
        }
      }

    }

  }

  static std::thread reporter;

  void addMemUsage(size_t bytes) {
    currentMemUsage += bytes;
    if(currentMemUsage > peakMemUsage) {
      peakMemUsage = currentMemUsage;
    }

    // now launch the periodic reporter
    // if(!reporter.joinable()) {
    //   reporter = std::thread(memUsagePeriodicReportInThread);
    // }
  }

  void removeMemUsage(size_t bytes) {
    currentMemUsage -= bytes;
    if(currentMemUsage == 0 && reporter.joinable()) {
      reporter.join();
    }

  }

  //! Class to manage and evaluate a field defined on a single grid.
  template<typename DataType, typename CoordinateType=tools::datatypes::strip_complex<DataType>>
  class Field : public std::enable_shared_from_this<Field<DataType, CoordinateType>> {
  public:
    using TGrid = const grids::Grid<CoordinateType>;
    using TPtrGrid = std::shared_ptr<TGrid>;
    using TData = std::vector<DataType>;
    using value_type = DataType;
    using ComplexType = tools::datatypes::ensure_complex<DataType>;

    using FourierManager = tools::numerics::fourier::FieldFourierManager<DataType, CoordinateType>;
    enum {
      x, y, z
    } DirectionType;


  protected:
    const TPtrGrid pGrid; //!< Pointer to the grid on which the field is defined.
    std::shared_ptr<FourierManager> fourierManager; //!< Class to handle Fourier transforms of this field.
    TData data; //!< Vector which stores the underlying data associated to the field
    bool fourier; //!< If true, then the field is regarded as being in Fourier space. Switched by Fourier transforms.

  public:
    //! Move constructor
    Field(Field<DataType, CoordinateType> &&move) : pGrid(move.pGrid),
                                                    fourier(move.fourier) {
      fourierManager = std::make_shared<FourierManager>(*this);
      std::swap(data, move.data);
      assert(data.size() == this->fourierManager->getRequiredDataSize());
    }

    //! Move operator
    auto & operator=(Field<DataType, CoordinateType> &&move) {
      assert(move.pGrid == pGrid);
      std::swap(data, move.data);
      assert(data.size() == this->fourierManager->getRequiredDataSize());
      fourier = move.fourier;
      return *this;
    }

    //! Copy constructor
    Field(const Field<DataType, CoordinateType> &copy)
      : std::enable_shared_from_this<Field<DataType, CoordinateType>>(),
        pGrid(copy.pGrid), data(copy.data),
        fourier(copy.fourier) {
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());
      addMemUsage(data.size() * sizeof(DataType));
    }

    //! Construct a field on the specified grid by moving the given data
    Field(TGrid &grid, TData &&dataVector, bool fourier = true) : pGrid(grid.shared_from_this()),
                                                                  data(std::move(dataVector)), fourier(fourier) {

      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());
      addMemUsage(data.size() * sizeof(DataType));

    }



    //! Construct a field on the specified grid by copying the given data
    Field(TGrid &grid, const TData &dataVector, bool fourier = true) : pGrid(grid.shared_from_this()),
                                                                       data(dataVector), fourier(fourier) {

      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());
      addMemUsage(data.size() * sizeof(DataType));

    }

    //! Construct a zero-filled field on the specified grid
    Field(TGrid &grid, bool fourier = true) : pGrid(grid.shared_from_this()),
                                              fourierManager(std::make_shared<FourierManager>(*this)),
                                              data(fourierManager->getRequiredDataSize(), 0),
                                              fourier(fourier) {
      addMemUsage(data.size() * sizeof(DataType));
    }

    virtual ~Field() {
      removeMemUsage(data.size() * sizeof(DataType));
    }

  public:

    //! \brief Returns a reference to the underlying grid
    TGrid &getGrid() const {
      return const_cast<TGrid &>(*pGrid);
    }

    //! Returns a reference to the data vector that stores the field
    TData &getDataVector() {
      return data;
    }

    //! Returns a constant reference to the data vector that stores the field.
    const TData &getDataVector() const {
      return data;
    }

    //! Returns a reference to the data vector storing the field.
    operator std::vector<DataType> &() {
      return getDataVector();
    }

    //! Returns a constant reference to the data vector storing the field.
    operator const std::vector<DataType> &() const {
      return getDataVector();
    }

    //! Evaluates the field at the grid point nearest to the supplied coordinate.
    DataType evaluateNearest(const Coordinate<CoordinateType> &location) const {
      auto offsetLower = pGrid->offsetLower;
      int x_p_0, y_p_0, z_p_0;


      // grid coordinates of parent cell that we're in
      x_p_0 = (int) floor(((location.x - offsetLower.x) / pGrid->cellSize));
      y_p_0 = (int) floor(((location.y - offsetLower.y) / pGrid->cellSize));
      z_p_0 = (int) floor(((location.z - offsetLower.z) / pGrid->cellSize));

      return (*this)[pGrid->getIndexFromCoordinate(Coordinate<int>(x_p_0, y_p_0, z_p_0))];
    }

    //! Multiply the field in-place by the provided field
    template<typename OtherDataType>
    void operator*=(const Field<OtherDataType, CoordinateType> & other) {
      // This operation only really make sense in real space
      assert(!other.isFourier() && !this->isFourier());
      // We need to use size3 rather than data.size() as `other` may not
      // be padded by FFTW (and hence it may be shorter).
      // To be on the safe side, we check that their grid have the same size.
      assert(other.getGrid().size3 == this->pGrid->size3);
      size_t N = this->pGrid->size3;
#pragma omp parallel for
      for(size_t i=0; i<N; i++) {
        data[i]*=other[i];
      }
    }

    //! Multiply the field in-place by the provided value
    template<typename OtherDataType>
    void operator*=(OtherDataType value) {
      size_t N = data.size();
#pragma omp parallel for
      for(size_t i=0; i<N; i++) {
        data[i]*=value;
      }
    }

    //! Add the provided field to this one in-place
    void operator+=(const Field<DataType, CoordinateType> & other) {
      size_t N = data.size();
#pragma omp parallel for
      for(size_t i=0; i<N; i++) {
        data[i]+=other[i];
      }
    }

    //! Subtract the provided field from this one in-place
    void operator-=(const Field<DataType, CoordinateType> & other) {
      size_t N = data.size();
#pragma omp parallel for
      for(size_t i=0; i<N; i++) {
        data[i]-=other[i];
      }
    }

    Field<DataType, CoordinateType> operator-(const Field<DataType, CoordinateType> & other) const {
      auto ret = Field<DataType, CoordinateType>(*this);
      ret-=other;
      return ret;
    }

    //! Add a multiple of the provided field to this one in-place
    void addScaled(const Field<DataType, CoordinateType> & other,
                   tools::datatypes::strip_complex<DataType> scale) {
      size_t N = data.size();
#pragma omp parallel for
      for(size_t i=0; i<N; i++) {
        data[i]+=scale*other[i];
      }
    }

    //! Single grid maximum. Only implemented for real space
    auto Maximum() const {
      assert(!isFourier());

      tools::datatypes::strip_complex<DataType> max=0;
      size_t N = data.size();
      
      for(size_t i=0; i<N; i++) {
        max = std::max(max, std::abs(data[i]));
      }
      return max;
    }

    //! Single grid inner product. Only implemented for real space.
    auto innerProduct(const Field<DataType, CoordinateType> & other) const {
      assert(!other.isFourier());
      assert(!isFourier());

      tools::datatypes::strip_complex<DataType> v=0;
      size_t N = data.size();

#pragma omp parallel for reduction(+:v)
      for(size_t i=0; i<N; i++) {
        v+=data[i]*other[i];
      }
      return v;
    }

    auto norm() const {
      auto norm2 = innerProduct(*this);
      return sqrt(norm2);
    }

    Field<DataType, CoordinateType> operator-() const {
      auto ret(*this);
      size_t N = data.size();
#pragma omp parallel for
      for(size_t i=0; i<N; i++) {
        ret[i]=-ret[i];
      }
      return ret;
    }

    void applyTransferFunction(const Field<DataType, CoordinateType> & covariance, double power) {
      using T = tools::datatypes::ensure_complex<DataType>;
      assert(this->isFourier());
      assert(covariance.isFourier());
      assert(&covariance.getGrid() == &this->getGrid());
      auto grid = this->getGrid();
      forEachFourierCellInt([&grid, this, &covariance, power]
                                    (T existingValue, int kx, int ky, int kz) {

        auto spec = covariance.getFourierCoefficient(kx,ky,kz).real();

        if(power!=1.0 && spec!=0.0)
          spec = pow(spec, power);

        T new_val = existingValue*spec;

        return new_val;
      });
    }



    //! Add the specified value to the field and the given location, using conjugate deinterpolation
    /*! For an explanation of what is meant by 'conjugate deinterpolation' see the
     * documentation for numerics::LocalUnitTricubicApproximation::getTransposeElementsForPosition
     */
    void deInterpolate(Coordinate<CoordinateType> location, DataType value) {

      // TODO: The below is not conjugate deinterpolation for the linear interpolation
      // implemented when CUBIC_INTERPOLATION is off. However, switching off
      // CUBIC_INTERPOLATION is not recommended, so fixing this inconsistency is low
      // priority.

      int x_p_0, y_p_0, z_p_0;
      CoordinateType dx, dy, dz;

      location -= pGrid->offsetLower;
      location = pGrid->wrapPoint(location);


      // grid coordinates of parent cell whose *centroid* (not corner) is to the bottom-left of our current point
      std::tie(x_p_0, y_p_0, z_p_0) = floor(location / pGrid->cellSize - 0.5);

      DataType valsForInterpolation[4][4][4];
      std::tie(dx,dy,dz) = (location / pGrid->cellSize - 0.5);
      dx-=x_p_0;
      dy-=y_p_0;
      dz-=z_p_0;

      // Caching would speed this up, since we anticipate repeated calls with the same dx,dy,dz (to numerical accuracy)
      numerics::LocalUnitTricubicApproximation<CoordinateType>::getTransposeElementsForPosition(dx,dy,dz,valsForInterpolation);



      // Figure out logistics of wrapping etc for indexing - more efficient than 4^3 calls to getIndexFromCoordinate
      int index_x[4], index_y[4], index_z[4];
      size_t gridSize = pGrid->size;
      size_t gridSize2 = gridSize*gridSize;
      auto getWrappedCoords = [gridSize](int coords[4], int key_cell_coord) {
        for(int i=0; i<4; ++i) {
          coords[i] = key_cell_coord + i - 1;
          if(coords[i]<0) coords[i]+=gridSize;
          if(coords[i]>gridSize) coords[i]-=gridSize;
        }
      };

      getWrappedCoords(index_x, x_p_0);
      getWrappedCoords(index_y, y_p_0);
      getWrappedCoords(index_z, z_p_0);

      size_t key_cell_index = pGrid->getIndexFromCoordinate({index_x[0], index_y[0], index_z[0]});

      auto & data = this->getDataVector();
      for(int i=0; i<4; ++i) {
        for(int j=0; j<4; ++j) {
          for(int k=0; k<4; ++k) {
            size_t this_index = key_cell_index +
                    (index_x[i] - index_x[0])*gridSize2 +
                    (index_y[j] - index_y[0])*gridSize +
                    (index_z[k] - index_z[0]);
            // above is an optimization of this_index = pGrid->getIndexFromCoordinate({x_p_0 - 1 + i, ...} )
            data[this_index] += value * valsForInterpolation[i][j][k];
          }
        }
      }

    }

    //! Evaluates the field at the specified co-ordinate using interpolation.
    DataType evaluateInterpolated(Coordinate<CoordinateType> location) const {
      int x_p_0, y_p_0, z_p_0;

      location -= pGrid->offsetLower;
      location = pGrid->wrapPoint(location);

      // grid coordinates of parent cell whose *centroid* (not corner) is to the bottom-left of our current point
      std::tie(x_p_0, y_p_0, z_p_0) = floor(location / pGrid->cellSize - 0.5);

#ifdef CUBIC_INTERPOLATION
      CoordinateType dx,dy,dz;

      // Work out the fractional displacement of our target point between the centroid of the cell identified above
      // and the next one along. dx, dy, dz will be between zero and one unless something goes badly wrong!
      std::tie(dx,dy,dz) = (location / pGrid->cellSize - 0.5);
      dx-=x_p_0;
      dy-=y_p_0;
      dz-=z_p_0;

      if(cache::enabled) {
        auto interp = getTricubicInterpolatorCached(x_p_0, y_p_0, z_p_0);
        return interp(dx, dy, dz);
      } else {
        auto interp = makeTricubicInterpolator(x_p_0, y_p_0, z_p_0);
        return interp(dx, dy, dz);
      }

#else

      int x_p_1, y_p_1, z_p_1;

      bool allowWrap = pGrid->coversFullSimulation();
      auto offsetLower = pGrid->offsetLower;

      // grid coordinates of top-right
      x_p_1 = x_p_0 + 1;
      y_p_1 = y_p_0 + 1;
      z_p_1 = z_p_0 + 1;

      // weights, which are the distance to the centre point of the
      // upper-right cell, in grid units (-> maximum 1)

      CoordinateType xw0, yw0, zw0, xw1, yw1, zw1;
      xw0 = ((CoordinateType) x_p_1 + 0.5) - (location.x / pGrid->cellSize);
      yw0 = ((CoordinateType) y_p_1 + 0.5) - (location.y / pGrid->cellSize);
      zw0 = ((CoordinateType) z_p_1 + 0.5) - (location.z / pGrid->cellSize);

      xw1 = 1. - xw0;
      yw1 = 1. - yw0;
      zw1 = 1. - zw0;

      assert(xw0 <= 1.0 && xw0 >= 0.0);

      int size_i = static_cast<int>(pGrid->size);

      if (allowWrap) {
        std::tie(x_p_0, y_p_0, z_p_0) = pGrid->wrapCoordinate({x_p_0, y_p_0, z_p_0});
        std::tie(x_p_1, y_p_1, z_p_1) = pGrid->wrapCoordinate({x_p_1, y_p_1, z_p_1});

      } else {
        // allow things on the boundary to 'saturate' value, but beyond boundary
        // is not acceptable
        assert(x_p_1 <= size_i);
        if (x_p_1 == size_i) x_p_1 = size_i - 1;
        assert(y_p_1 <= size_i);
        if (y_p_1 == size_i) y_p_1 = size_i - 1;
        assert(z_p_1 <= size_i);
        if (z_p_1 == size_i) z_p_1 = size_i - 1;
        assert(x_p_0 >= -1);
        if (x_p_0 == -1) x_p_0 = 0;
        assert(y_p_0 >= -1);
        if (y_p_0 == -1) y_p_0 = 0;
        assert(z_p_0 >= -1);
        if (z_p_0 == -1) z_p_0 = 0;

      }

      // at this point every evaluation about to take place should be at a proper grid location
      assert(x_p_0 < size_i && x_p_0 >= 0 && x_p_1 < size_i && x_p_1 >= 0);
      assert(y_p_0 < size_i && y_p_0 >= 0 && y_p_1 < size_i && y_p_1 >= 0);
      assert(z_p_0 < size_i && z_p_0 >= 0 && z_p_1 < size_i && z_p_1 >= 0);


      return xw0 * yw0 * zw1 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_0, y_p_0, z_p_1)] +
             xw1 * yw0 * zw1 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_1, y_p_0, z_p_1)] +
             xw0 * yw1 * zw1 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_0, y_p_1, z_p_1)] +
             xw1 * yw1 * zw1 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_1, y_p_1, z_p_1)] +
             xw0 * yw0 * zw0 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_0, y_p_0, z_p_0)] +
             xw1 * yw0 * zw0 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_1, y_p_0, z_p_0)] +
             xw0 * yw1 * zw0 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_0, y_p_1, z_p_0)] +
             xw1 * yw1 * zw0 * (*this)[pGrid->getIndexFromCoordinateNoWrap(x_p_1, y_p_1, z_p_0)];

#endif

    }



  protected:


    const numerics::LocalUnitTricubicApproximation<DataType> getTricubicInterpolatorCached(int x_p_0, int y_p_0, int z_p_0) const {
      assert(cache::enabled);
      auto key = std::make_tuple(x_p_0, y_p_0, z_p_0, static_cast<const void *>(this));
      auto result = cache::cachedInterpolators<DataType>->get(key);
      if (result == boost::none) {
        cache::cacheMisses += 1;
        cache::cachedInterpolators<DataType>->insert(key, makeTricubicInterpolator(x_p_0, y_p_0, z_p_0));
        return cache::cachedInterpolators<DataType>->get(key).get();
      } else {
        cache::cacheHits += 1;
        return result.get();
      }
    }

    numerics::LocalUnitTricubicApproximation<DataType> makeTricubicInterpolator(int x_p_0, int y_p_0, int z_p_0) const {
      assert(!this->isFourier());
      DataType valsForInterpolation[4][4][4];
      for(int i=-1; i<3; ++i) {
        for(int j=-1; j<3; ++j) {
          for(int k=-1; k<3; ++k) {
            Coordinate<int> coord {x_p_0+i, y_p_0+j, z_p_0+k};
            if(pGrid->size == pGrid->simEquivalentSize) {
              coord = pGrid->wrapCoordinate(coord);
            } else {
              // Repeat values at the boundary. Note this leads to artefacts in the interpolation. Warnings
              // about using results within a few cells of the boundary are issued to the user elsewhere.
              coord = pGrid->clampCoordinate(coord);
            }
            valsForInterpolation[i+1][j+1][k+1] = (*this)[pGrid->getIndexFromCoordinateNoWrap(coord)];
          }
        }
      }
      return numerics::LocalUnitTricubicApproximation<DataType>(valsForInterpolation);
    }

  public:

    //! Returns a constant reference to the value of the field at grid index i
    const DataType &operator[](size_t i) const {
      return data[i];
    }

    //! Returns a reference to the value of the field at grid index i
    DataType &operator[](size_t i) {
      return data[i];
    }

    /*
     * FOURIER TRANSFORM SUPPORT
     */

    //! Converts the field to the same co-ordinates (Fourier or Real space) as the specified field.
    template<typename T>
    void matchFourier(const T &other) {
      if (other.isFourier())
        toFourier();
      else
        toReal();
    }

    //! Returns true if the field is in Fourier space.
    bool isFourier() const {
      return fourier;
    }

    //! Asserts whether the field is in Fourier space, without actually applying any transform
    void setFourier(bool fourier) {
      this->fourier = fourier;
    }

    //! \brief Converts the field to Fourier space
    /*!
        Does nothing if already in Fourier space.
    */
    void toFourier() {
      if (fourier) return;
      fourierManager->performTransform();
      assert(fourier);
    }

    //! \brief Converts the field to real space
    /*!
        Does nothing if already in real space.
    */
    void toReal() {
      if (!fourier) return;
      fourierManager->performTransform();
      assert(!fourier);
    }

    //! Returns the value of the field in Fourier space at the specified Fourier mode
    //! For efficiency, does not check whether the field is actually stored in Fourier space first.
    ComplexType getFourierCoefficient(int kx, int ky, int kz) const {
      return fourierManager->getFourierCoefficient(kx, ky, kz);
    }

    //! Generate a set of three Fourier fields from a function of k and the Fourier space field, supplied as an argument.
    template<typename... Args>
    auto generateNewFourierFields(Args &&... args) {
      return fourierManager->generateNewFourierFields(args...);
    }

    //! Iterate (potentially in parallel) over each Fourier cell.
    /*!
     * The passed function takes arguments (value, kx, ky, kz) where value is the Fourier coeff value
     * at k-mode kx, ky, kz, and kx,ky,kz are the modes in comoving (h/Mpc) coordinates.
     * If the function returns a value, the Fourier coeff in that cell is updated accordingly.
     */
    template<typename... Args>
    void forEachFourierCell(Args &&... args) {
      fourierManager->forEachFourierCell(args...);
    }

    //! Overload of forEachFourierCell that cannot modify the field
    template<typename... Args>
    void forEachFourierCell(Args &&... args) const {
      fourierManager->forEachFourierCell(args...);
    }

    //! Iterate over Fourier cells with a function taking integer kx,ky,kz arguments.
    template<typename... Args>
    void forEachFourierCellInt(Args &&... args) const {
      fourierManager->forEachFourierCellInt(args...);
    }

    //! Iterate (potentially in parallel) over each Fourier cell, accumulate a complex number over each Fourier cell
    /*!
    The passed function takes arguments (value, kx, ky, kz). The return value is accumulated.
    */
    template<typename... Args>
    auto accumulateForEachFourierCell(Args &&... args) const {
      return fourierManager->accumulateForEachFourierCell(args...);
    }

    //! \brief Sets the value of the field in Fourier space, at the specified mode.
    /*!
    \param kx - integer kx mode
    \param ky - integer ky mode
    \param kz - integer kz mode
    \param value - value to set the field to at this mode
    */
    void setFourierCoefficient(int kx, int ky, int kz, const ComplexType &value) {
      fourierManager->setFourierCoefficient(kx, ky, kz, value);
    }

    /*! Ensure that the Fourier modes in this field are correctly 'mirrored', i.e. where there is duplication,
     * consistent values are stored. This likely only becomes an issue for real FFTs.
     */
    void ensureFourierModesAreMirrored() const {
      assert(isFourier());

      // Logically this is a const operation, even though actually it will potentially change internal state. So our
      // externally visible method is const, but the actual manipulation is non-const.
      const_cast<FourierManager &>(*fourierManager).ensureFourierModesAreMirrored();
    }

    //! Apply a Fourier space filter that suppresses the field at some k
    void applyFilter(const filters::Filter<CoordinateType> &filter) {
      forEachFourierCell([&filter](ComplexType current_value, CoordinateType kx, CoordinateType ky, CoordinateType kz) {
        CoordinateType k = sqrt(double(kx * kx + ky * ky + kz * kz));
        return current_value * filter(k);
      });
    }

    /*! \brief Apply a Fourier space filter only in the specified window
     *
     * If windowFirst is true, the filter is applied AFTER the masking takes place; otherwise the order is reversed.
     */
    void applyFilterInWindow(const filters::Filter<CoordinateType> &filter, const Window<CoordinateType> & window,
                             bool windowFirst) {
      using tools::numerics::operator+=;

      auto inWindow = this->copy();

      if(windowFirst) {
        inWindow->setZeroOutsideWindow(window);
        inWindow->toFourier();
        inWindow->applyFilter(filter);
      } else {
        inWindow->toFourier();
        inWindow->applyFilter(filter);
        inWindow->setZeroOutsideWindow(window);
      }

      this->setZeroInsideWindow(window);
      if(inWindow->isFourier())
        this->toFourier();
      else
        this->toReal();
      this->data += inWindow->data;

    }

    void setZeroInsideWindow(const Window<CoordinateType> & window) {
      toReal();

#pragma omp parallel for default(none) shared(window)
      for(size_t i=0; i<this->pGrid->size3; ++i) {
        if(window.contains(this->pGrid->getCentroidFromIndex(i)))
          (*this)[i]=0;
      }
    }

    void setZeroOutsideWindow(const Window<CoordinateType> & window) {
      toReal();

#pragma omp parallel for default(none) shared(window)
      for(size_t i=0; i<this->pGrid->size3; ++i) {
        if(!window.contains(this->pGrid->getCentroidFromIndex(i)))
          (*this)[i]=0;
      }
    }

    //! Adds the supplied field to this one, even if it is defined using a different grid.
    /*!
     * Requires the source field to be in real (rather than Fourier) space.
     *
     *  The addition performs interpolation if necessary. It does not modify the grid
     *  that is being added to this one.
    */
    void addFieldFromDifferentGrid(const Field<DataType, CoordinateType> &source) {
      assert(!source.isFourier());
      toReal();
      TPtrGrid pSourceProxyGrid = source.getGrid().makeProxyGridToMatch(getGrid());

      auto evaluator = makeEvaluator(source, *pSourceProxyGrid);
      evaluator->addTo(*this);

    }


    //! Add a field defined on a different grid to this one, applying the specified filter while adding
    void addFieldFromDifferentGridWithFilter(const Field<DataType, CoordinateType> & source,
                                             const filters::Filter<CoordinateType> & filter) {

#ifdef FILTER_ON_COARSE_GRID
      // Version for compatibility with pre-Dec 2016 output
      // Applies filter BEFORE interpolating onto the fine grid, which results in more pixel window function
      // artefacts but (presumably?) fewer artefacts from the grid-level window function

      auto temporaryField = std::make_shared<Field<DataType, CoordinateType>>(source);
      // counterintuitively requires a shared_ptr due to innards of addFieldFromDifferentGrid using shared_from_this

      temporaryField->applyFilter(filter);
      temporaryField->toReal();
      this->toReal();
      this->addFieldFromDifferentGrid(*temporaryField);

#else
      assert(!source.isFourier());

      auto temporaryField = std::make_shared<Field<DataType, CoordinateType>>(getGrid(), false);

      temporaryField->addFieldFromDifferentGrid(source);
      temporaryField->applyFilter(filter);

      this->matchFourier(*temporaryField); // expect that the temporary field is now stored in Fourier space

      const auto & temporaryFieldData = temporaryField->getDataVector();

#pragma omp parallel for schedule(static) default(none)  shared(data, temporaryFieldData)
      for (size_t i = 0; i < temporaryFieldData.size(); ++i) {
        data[i] += temporaryFieldData[i];
      }
#endif
    } 

    //! addFieldFromDifferentGrid, automatically converting the source field from Fourier space if required
    void addFieldFromDifferentGrid(Field<DataType, CoordinateType> &source) {
#ifndef FILTER_ON_COARSE_GRID
      source.toReal();
#endif
      addFieldFromDifferentGrid(const_cast<const Field<DataType, CoordinateType> &>(source));
    }

      //! addFieldFromDifferentGridWithFilter, automatically converting the source field from Fourier space if required
      void addFieldFromDifferentGridWithFilter(Field<DataType, CoordinateType> &source,
                                               const filters::Filter<CoordinateType> & filter) {
#ifndef FILTER_ON_COARSE_GRID
        source.toReal();
#endif
        addFieldFromDifferentGridWithFilter(const_cast<const Field<DataType, CoordinateType> &>(source), filter);
      }

    //! Outputs the field as a numpy array to the specified filename.
    void dumpGridData(std::string filename) const {
      int n = static_cast<int>(getGrid().size);
      const int dim[3] = {n, n, n};
      io::numpy::SaveArrayAsNumpy(filename, false, 3, dim, data.data());
    }

    //! Loads field data from the specified numpy file, if this is possible.
    void loadGridData(std::string filename) {
      int n = static_cast<int>(getGrid().size);
      int n0, n1, n2;
      io::numpy::LoadArrayFromNumpy(filename, n0, n1, n2, data);
      if (n0 != n || n1 != n || n2 != n) {
        throw std::runtime_error("Incorrect size for imported numpy array");
      }
      assert(data.size() == getGrid().size3);
      data.resize(fourierManager->getRequiredDataSize());
    }

    auto copy() const {
      return std::make_shared<Field<DataType,CoordinateType>>(*this);
    }


  };

  //! Converts field from one storage type to another, likely from complex -> real internal representation.
  //! Note that convertField is only retained for future debugging.
  template<typename TargetDataType, typename SourceDataType, typename CoordinateType>
  std::shared_ptr<Field<TargetDataType, CoordinateType>>
  convertField(Field<SourceDataType, CoordinateType> &field) {
    assert(!field.isFourier());
    auto newField = make_shared<Field<TargetDataType, CoordinateType>>(field.getGrid(), field.isFourier());
    auto &newData = newField->getDataVector();
    auto &originalData = field.getDataVector();
    size_t size = field.getGrid().size3;

#pragma omp parallel for schedule(static) default(none) shared(newData, originalData)
    for (size_t i = 0; i < size; ++i) {
      newData[i] = tools::datatypes::real_part_if_complex(originalData[i]);
    }

    return newField;
  };

  //! Specialisation of convertField for the case where no conversion is necessary.
  //! Note that convertField is only retained for future debugging.
  template<typename TargetDataType, typename CoordinateType>
  std::shared_ptr<Field<TargetDataType, CoordinateType>>
  convertField(Field<TargetDataType, CoordinateType> &field) {
    return field.shared_from_this();
  };


}


#endif
