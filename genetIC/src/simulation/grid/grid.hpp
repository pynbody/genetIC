#ifndef __GRID_HPP
#define __GRID_HPP

#include <cassert>
#include <set>
#include <type_traits>
#include <memory>
#include <vector>
#include <complex>
#include "src/tools/numerics/fourier.hpp"
#include "src/simulation/coordinate.hpp"
#include "src/tools/progress/progress.hpp"
#include "src/tools/util_functions.hpp"
#include "src/tools/data_types/complex.hpp"
#include "src/simulation/grid/virtualgrid.hpp"

using std::complex;
using std::vector;
using std::make_shared;

/*!
    \namespace grids
    \brief Define a grid with given size and resolution. This a building block for fields and constraints.

    Grids contain cells which can be referred to by either

     * a coordinate, which is defined as integers with (0,0,0) being the bottom-left corner, (size-1,size-1,size-1)
       being the top right;
     * an index, which is defined as size_t from 0 to size^3
     * a point/centroid, which is defined as a floating point (type T) triple. The bottom-left corner of the grid is given by
       offsetLower, and the top right by offsetLower + thisGridSize.

 */

namespace grids {

  template<typename T>
  class Grid : public std::enable_shared_from_this<Grid<T>> {
  public:

    using GridPtrType = std::shared_ptr<Grid<T>>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;

  private:
    T kMin; /*!< Fundamental mode of the box */
    std::vector<size_t> flags;  /*!< Flagged cells on this grid */

  public:
    const T periodicDomainSize;
    const T thisGridSize;   /*!< Grid (one side) size in Mpc */
    const T cellSize;
    const Coordinate<T> offsetLower; ///< Coordinate of the pixel edge of lower corner of the grid
    const size_t size; ///<the number of cells on a side
    const size_t size2; ///< the number of cells on a face
    const size_t size3; ///< the total number of cells in the grid cube
    const size_t simEquivalentSize; ///< the total number of cells on the box side of a simulation if it were all at this resolution
    const T cellMassFrac; ///< the fraction of mass of the full simulation in a single cell of this grid
    const T cellSofteningScale; ///< normally 1.0; scales softening relative to dx

    Grid(T simsize, size_t n, T dx = 1.0, T x0 = 0.0, T y0 = 0.0, T z0 = 0.0,
         T massFrac = 0.0, T softScale = 1.0) :
        periodicDomainSize(simsize), thisGridSize(dx * n),
        cellSize(dx), offsetLower(x0, y0, z0),
        size(n), size2(n * n), size3(n * n * n),
        simEquivalentSize((unsigned) tools::getRatioAndAssertInteger(simsize, dx)),
        cellMassFrac(massFrac == 0.0 ? pow(dx / simsize, 3.0) : massFrac),
        cellSofteningScale(softScale) {
      setKmin();
    }


    explicit Grid(size_t n) : periodicDomainSize(0), thisGridSize(n),
                     cellSize(1.0), offsetLower(0, 0, 0),
                     size(n), size2(n * n), size3(n * n * n), simEquivalentSize(0), cellMassFrac(0.0),
                     cellSofteningScale(1.0) {
      setKmin();
    }


  protected:
    void setKmin() {
      kMin = 2. * M_PI / thisGridSize;
    }

  public:

    bool coversFullSimulation() const {
      return periodicDomainSize == thisGridSize;
    }

    T getFourierKmin() const {
      return kMin;
    }

    virtual void debugInfo(std::ostream &s) const {
      s << "Grid of side " << size << " address " << this << "; " << this->flags.size() << " cells marked";
    }

    int getEffectiveSimulationSize() const {
      return tools::getRatioAndAssertInteger(periodicDomainSize, cellSize);
    }

    /*****************************
    * Methods dealing with the creation of virtual grids and relationships between grids.
    ******************************/

    bool pointsToAnyGrid(std::vector<std::shared_ptr<Grid<T>>> grids) {
      for (auto g: grids) {
        if (pointsToGrid(g.get()))
          return true;
      }
      return false;
    }

    virtual bool isUpsampledOrDownsampled(){
      return false;
    }


    /*! Return true if this grid has a known relationship to pOther, in the sense that a field defined on pOther
    * could be evaluated on this grid. */
    virtual bool pointsToGrid(const Grid *pOther) const {
      return this == pOther;
    }


    //! Make a VirtualGrid pointing to this grid, but with a resolution and range matching target.
    GridPtrType makeProxyGridToMatch(const Grid<T> &target) const {
      GridPtrType proxy = std::const_pointer_cast<Grid<T>>(this->shared_from_this());
      if (target.cellSize > cellSize) {
        size_t ratio = tools::getRatioAndAssertPositiveInteger(target.cellSize, cellSize);
        proxy = std::make_shared<SubSampleGrid<T>>(proxy, ratio);
      } else if (target.cellSize < cellSize) {
        size_t ratio = tools::getRatioAndAssertPositiveInteger(cellSize, target.cellSize);
        proxy = std::make_shared<SuperSampleGrid<T>>(proxy, ratio);
      }

      if (target.offsetLower != offsetLower || target.size != proxy->size) {
        auto relativeOffset = target.offsetLower - offsetLower;
        proxy = std::make_shared<SectionOfGrid<T>>(proxy,
                                                   tools::getRatioAndAssertInteger(relativeOffset.x,
                                                                                   proxy->cellSize),
                                                   tools::getRatioAndAssertInteger(relativeOffset.y,
                                                                                   proxy->cellSize),
                                                   tools::getRatioAndAssertInteger(relativeOffset.z,
                                                                                   proxy->cellSize),
                                                   target.size);
      }


      return proxy;
    }

    virtual std::shared_ptr<Grid<T>> makeScaledMassVersion(T massRatio) {
      return std::make_shared<MassScaledGrid<T>>(this->shared_from_this(), massRatio);
    }

    /*****************************
     * Methods dealing with flagging cells on the grid
     ******************************/

    virtual void getFlaggedCells(std::vector<size_t> &targetArray) const {
      targetArray.insert(targetArray.end(), flags.begin(), flags.end());
    }

    virtual void flagCells(const std::vector<size_t> &sourceArray) {
      std::vector<size_t> newFlags(flags.size() + sourceArray.size());
      auto end = std::set_union(flags.begin(), flags.end(), sourceArray.begin(), sourceArray.end(), newFlags.begin());
      newFlags.resize(end - newFlags.begin());
      flags = std::move(newFlags);
      tools::sortAndEraseDuplicate(flags);
    }

    virtual void expandFlaggedRegionInDirection(const Coordinate<int> &step) {
      size_t old_size = flags.size();
      flags.resize(old_size*3);
      for (size_t i=0; i<old_size; ++i) {
        size_t original_cell_id = flags[i];
        flags[i+old_size] = this->getIndexFromIndexAndStep(original_cell_id, step);
        flags[i+old_size*2] = this->getIndexFromIndexAndStep(original_cell_id, -step);
      }
      tools::sortAndEraseDuplicate(flags);
    }

    virtual void expandFlaggedRegion(size_t ncells=1) {
      for(size_t i=0; i<ncells; i++) {
        expandFlaggedRegionInDirection({0,0,1});
        expandFlaggedRegionInDirection({0,1,0});
        expandFlaggedRegionInDirection({1,0,0});
      }
    }

    virtual void unflagAllCells() {
      flags.clear();
    }

    virtual size_t numFlaggedCells() const {
      return flags.size();
    }

    Coordinate<T> getFlaggedCellsCentre(){
      return this->getCentreWrapped(this->flags);
    }

    int getFlaggedCellsSize(){
      if (this->numFlaggedCells() > 0 ) {
          Window<int> flaggedWindow(this->getEffectiveSimulationSize(),
                                    this->getCoordinateFromIndex(this->flags[0]));
          for (auto cell_id : this->flags) {
              flaggedWindow.expandToInclude(this->getCoordinateFromIndex(cell_id));
          }
          return flaggedWindow.getMaximumDimension();
      }
      else {
              throw std::runtime_error("Trying to obtain physical size of an empty flagged region");
          }
    }

    T getFlaggedCellsPhysicalSize(){
        return T(this->getFlaggedCellsSize()) * this->cellSize;
    };

  protected:

    static void upscaleCellFlagVector(const std::vector<size_t> sourceArray,
                                      std::vector<size_t> &targetArray,
                                      const Grid<T> *source,
                                      const Grid<T> *target) {

      assert(target->size >= source->size);
      assert((target->size) % (source->size) == 0);
      int factor = int(target->size / source->size);
      targetArray.clear();

      for (auto id: sourceArray) {
        auto coord = source->getCoordinateFromIndex(id);
        iterateOverCube<int>(
            coord * factor, coord * factor + factor,
            [&targetArray, &target](const Coordinate<int> &subCoord) {
              targetArray.push_back(target->getIndexFromCoordinateNoWrap(subCoord));
            }
        );
      }
    }

    static void downscaleCellFlagVector(const std::vector<size_t> sourceArray,
                                        std::vector<size_t> &targetArray,
                                        const Grid<T> *source,
                                        const Grid<T> *target) {


      targetArray.clear();
      targetArray.resize(sourceArray.size());

      assert(source->size >= target->size);
      assert(source->offsetLower == target->offsetLower);

      size_t factor = tools::getRatioAndAssertPositiveInteger(target->cellSize, source->cellSize);

      // it's not clear that the following parallelisation actually speeds much up
#pragma omp parallel for
      for (size_t i = 0; i < sourceArray.size(); ++i) {
        size_t id = sourceArray[i];
        auto coord = source->getCoordinateFromIndex(id);
        targetArray[i] = target->getIndexFromCoordinateNoWrap(coord / factor);
      }

      // this sort seems to be the slowest step. In C++17 we can make it parallel... or is there a
      // better overall algorithm?
      tools::sortAndEraseDuplicate(targetArray);
    }


    /*****************************
    * Methods dealing with all sorts of coordinate calculations
    ******************************/

  public:
    T getWrappedOffset(T x0, T x1) const {
      T result = x0 - x1;
      if (result > periodicDomainSize / 2) {
        result -= periodicDomainSize;
      }
      if (result < -periodicDomainSize / 2) {
        result += periodicDomainSize;
      }
      return result;
    }

    Coordinate<T> getWrappedOffset(const Coordinate<T> a, const Coordinate<T> b) const {
      return Coordinate<T>(getWrappedOffset(a.x, b.x), getWrappedOffset(a.y, b.y), getWrappedOffset(a.z, b.z));
    }

    //! Calculate the centre in box coordinate of a vector of ids
    /*! The underlying assumption of this method is that the centering is done on the coarse grid.
     * Centering on zoom grids is not taken care off.
     */
    Coordinate<T> const getCentreWrapped(const std::vector<size_t>& vector_ids){
      if(vector_ids.empty()){
        throw std::runtime_error("Cannot calculate the center of an empty region");
      }

      T runningx = 0.0;
      T runningy = 0.0;
      T runningz = 0.0;

      auto p0_location = this->getCentroidFromIndex(vector_ids[0]);

      // Calculate the wrapped mean wrto to cell 0
      for (size_t i = 1; i <vector_ids.size(); i++) {
        size_t id = vector_ids[i];
        auto pi_location = this->getCentroidFromIndex(id);
        runningx += this->getWrappedOffset(pi_location.x, p0_location.x);
        runningy += this->getWrappedOffset(pi_location.y, p0_location.y);
        runningz += this->getWrappedOffset(pi_location.z, p0_location.z);
      }
      runningx /= vector_ids.size();
      runningy /= vector_ids.size();
      runningz /= vector_ids.size();

      // Add back cell 0 and wrap if needed
      runningx += p0_location.x;
      runningy += p0_location.y;
      runningz += p0_location.z;
      return this->wrapPoint(Coordinate<T>(runningx, runningy, runningz));
    }

    //! True if point in physical coordinates is on this grid
    virtual bool containsPoint(const Coordinate<T> &coord) const {
      return Window<T>(periodicDomainSize, offsetLower, offsetLower + thisGridSize).contains(coord);
    }

    //! True if point in physical coordinates is on this grid and not too close to the border
    /*!
     * @param safety Exclude "safety" number of pixels at the edge of the box
     */
    virtual bool containsPointWithBorderSafety(const Coordinate<T> &coord, int safety ) const{
      if(safety < 1){
        throw std::runtime_error("Safety number of pixels must be at least one");
      }

      return Window<T>(periodicDomainSize, offsetLower,
                       offsetLower + thisGridSize).containsWithBorderSafety(coord, safety * cellSize);
    }

    Coordinate<T> wrapPoint(Coordinate<T> pos) const {
      pos.x = fmod(pos.x, periodicDomainSize);
      if (pos.x < 0) pos.x += periodicDomainSize;
      pos.y = fmod(pos.y, periodicDomainSize);
      if (pos.y < 0) pos.y += periodicDomainSize;
      pos.z = fmod(pos.z, periodicDomainSize);
      if (pos.z < 0) pos.z += periodicDomainSize;
      return pos;
    }

    //! True if cell with pixel coordinates is on this grid
    /*! Does not take into account offset or physical coordinates
     */
    virtual bool containsCellWithCoordinate(Coordinate<int> coord) const {
      return coord.x >= 0 && coord.y >= 0 && coord.z >= 0 &&
             (unsigned) coord.x < size && (unsigned) coord.y < size && (unsigned) coord.z < size;
    }

    //! True if cell number is less than Ncell cubed
    virtual bool containsCell(size_t i) const {
      return i < size3;
    }


    size_t getIndexFromIndexAndStep(size_t index, const Coordinate<int> &step) const {
      auto coord = getCoordinateFromIndex(index);
      coord += step;
      return this->getIndexFromCoordinate(coord); // N.B. does wrapping inside getIndex
    }


    /*! Wrap the coordinate such that it lies within [0,size) of the base grid if this is possible.
     *
     * Note that for efficiency this routine only "corrects" coordinates within one boxsize of the fundamental domain.
     */
    Coordinate<int> wrapCoordinate(Coordinate<int> coord) const {
      if (coord.x > (signed) simEquivalentSize - 1) coord.x -= simEquivalentSize;
      if (coord.y > (signed) simEquivalentSize - 1) coord.y -= simEquivalentSize;
      if (coord.z > (signed) simEquivalentSize - 1) coord.z -= simEquivalentSize;
      if (coord.x < 0) coord.x += simEquivalentSize;
      if (coord.y < 0) coord.y += simEquivalentSize;
      if (coord.z < 0) coord.z += simEquivalentSize;
      return coord;
    }


    virtual size_t getIndexFromCoordinate(Coordinate<int> coord) const {
      coord = wrapCoordinate(coord);
      return getIndexFromCoordinateNoWrap(coord);
    }

    virtual size_t getIndexFromCoordinateNoWrap(size_t x, size_t y, size_t z) const {
      size_t index = (x * size + y) * size + z;
      assert(this->containsCell(index));
      return index;
    }

    virtual size_t getIndexFromCoordinateNoWrap(int x, int y, int z) const {

#ifdef SAFER_SLOWER
      if(x<0 || x>=size || y<0 || y>=size || z<0 || z>=size)
          throw std::runtime_error("Grid index out of range in getIndexNoWrap");
#endif
      auto index = size_t(x * size + y) * size + z;
      assert(this->containsCell(index));
      return index;
    }

     virtual size_t getIndexFromCoordinateNoWrap(const Coordinate<int> &coordinate) const {
      return getIndexFromCoordinateNoWrap(coordinate.x, coordinate.y, coordinate.z);
    }

    //! Returns cell id in pixel coordinates
    virtual Coordinate<int> getCoordinateFromIndex(size_t id) const {
      size_t x, y;

      if ((unsigned) id >= size3) {
        throw std::runtime_error("Index out of range");
      }

      // The following implementation is a little faster than using the
      // modulo operator.
      x = id / size2;
      id -= x * size2;
      y = id / size;
      id -= y * size;

      auto coord = Coordinate<int>(int(x), int(y), int(id));
      assert(this->containsCellWithCoordinate(coord));
      return coord;
    }

    //! Returns coordinate of centre of cell id, in physical box coordinates
    /*! Takes into account grid offsets wrt base grid, pixel size etc
     */
    virtual Coordinate<T> getCentroidFromIndex(size_t id) const {
      Coordinate<int> coord = getCoordinateFromIndex(id);
      return getCentroidFromCoordinate(coord);
    }


    virtual Coordinate<T> getCentroidFromCoordinate(const Coordinate<int> &coord) const {
      Coordinate<T> result(coord);
      result *= cellSize;
      result += offsetLower;
      result += cellSize / 2;
      assert(this->containsPoint(result));
      return result;
    }

    virtual size_t getIndexFromPoint(Coordinate<T> point) const {
//      auto coords = floor(wrapPoint(point - offsetLower - cellSize / 2) / cellSize);
       auto coords = floor(wrapPoint(point - offsetLower) / cellSize);
      assert(this->containsCellWithCoordinate(coords));
      return getIndexFromCoordinateNoWrap(coords);
    }

    virtual Coordinate<int> getCoordinateFromPoint(Coordinate<T> point) const {
      return this->getCoordinateFromIndex(this->getIndexFromPoint(point));
    }

    /*****************************
    * Methods dealing with insertion of new ids
    ******************************/

    void appendIdsInCubeToVector(T x0c, T y0c, T z0c, T dxc, vector<size_t> &ids) {
      size_t offset = ids.size();
      int added_size = std::round(dxc / cellSize);
      added_size *= added_size * added_size;
      ids.resize(offset + added_size);
      insertCubeIdsIntoVector(x0c, y0c, z0c, dxc, ids.begin() + offset);
    }

    void insertCubeIdsIntoVector(T x0c, T y0c, T z0c, T dxc, vector<size_t>::iterator start) {
      // return all the grid IDs whose centres lie within the specified cube

      std::tie(x0c, y0c, z0c) = wrapPoint(Coordinate<T>(x0c, y0c, z0c) - offsetLower);

      int xa = ((int) floor((x0c - dxc / 2 + cellSize / 2) / cellSize));
      int ya = ((int) floor((y0c - dxc / 2 + cellSize / 2) / cellSize));
      int za = ((int) floor((z0c - dxc / 2 + cellSize / 2) / cellSize));

      int xb = ((int) floor((x0c + dxc / 2 - cellSize / 2) / cellSize));
      int yb = ((int) floor((y0c + dxc / 2 - cellSize / 2) / cellSize));
      int zb = ((int) floor((z0c + dxc / 2 - cellSize / 2) / cellSize));

      iterateOverCube<int>(Coordinate<int>(xa, ya, za),
                           Coordinate<int>(xb, yb, zb) + 1,
                           [&start, this](const Coordinate<int> &cellCoord) {
                             (*start) = getIndexFromCoordinate(cellCoord);
                             assert(*start < size3);
                             ++start;
                           });

    }
  };
}

#endif

