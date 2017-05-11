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

 */

namespace grids {

  template<typename T>
  class Grid : public std::enable_shared_from_this<Grid<T>> {
  public:

    using GridPtrType = std::shared_ptr<Grid<T>>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;

  private:
    T kMin, kMinSquared;
    std::vector<size_t> flags;

  public:
    const T simsize, boxsize, dx;
    const Coordinate<T> offsetLower;
    const size_t size; ///<the number of cells on a side
    const size_t size2; ///< the number of cells on a face
    const size_t size3; ///< the total number of cells in the grid cube
    const size_t simEquivalentSize; ///< the total number of cells on the box side of a simulation if it were all at this resolution
    const T cellMassFrac; ///< the fraction of mass of the full simulation in a single cell of this grid
    const T cellSofteningScale; ///< normally 1.0; scales softening relative to dx

    Grid(T simsize, size_t n, T dx = 1.0, T x0 = 0.0, T y0 = 0.0, T z0 = 0.0,
         T massFrac = 0.0, T softScale = 1.0) :
      simsize(simsize), boxsize(dx * n),
      dx(dx), offsetLower(x0, y0, z0),
      size(n), size2(n * n), size3(n * n * n),
      simEquivalentSize((unsigned) tools::getRatioAndAssertInteger(simsize,dx)),
      cellMassFrac(massFrac == 0.0 ? pow(dx / simsize, 3.0) : massFrac),
      cellSofteningScale(softScale) {
      // cerr << "Grid ctor " << this <<  endl;
      setKmin();
    }


    Grid(size_t n) : simsize(0), boxsize(n),
                     dx(1.0), offsetLower(0, 0, 0),
                     size(n), size2(n * n), size3(n * n * n), simEquivalentSize(0), cellMassFrac(0.0),
                     cellSofteningScale(1.0) {
      // cerr << "Grid ctor size-only" << endl;
      setKmin();
    }

    virtual ~Grid() {
      // cerr << "~Grid " << this << endl;
    }

  protected:
    void setKmin() {
      kMin = 2. * M_PI / boxsize;
      kMinSquared = kMin * kMin;
    }

  public:

    bool coversFullSimulation() const {
      return simsize==boxsize;
    }

    T getWrappedDelta(T x0, T x1) const {
      T result = x0 - x1;
      if (result > simsize / 2) {
        result -= simsize;
      }
      if (result < -simsize / 2) {
        result += simsize;
      }
      return result;
    }

    Coordinate<T> getWrappedDelta(const Coordinate<T> a, const Coordinate<T> b) const {
      return Coordinate<T>(getWrappedDelta(a.x,b.x), getWrappedDelta(a.y, b.y), getWrappedDelta(a.z, b.z));
    }

    GridPtrType makeProxyGridToMatch(const Grid<T> &target) const {
      GridPtrType proxy = std::const_pointer_cast<Grid<T>>(this->shared_from_this());
      if (target.dx > dx) {
        size_t ratio = tools::getRatioAndAssertPositiveInteger(target.dx, dx);
        proxy = std::make_shared<SubSampleGrid<T>>(proxy, ratio);
      } else if (target.dx < dx) {
        size_t ratio = tools::getRatioAndAssertPositiveInteger(dx, target.dx);
        proxy = std::make_shared<SuperSampleGrid<T>>(proxy, ratio);
      }

      if (target.offsetLower != offsetLower || target.size != proxy->size) {
        auto relativeOffset = target.offsetLower-offsetLower;
        proxy = std::make_shared<SectionOfGrid<T>>(proxy,
                                                   tools::getRatioAndAssertInteger(relativeOffset.x,
                                                                                   proxy->dx),
                                                   tools::getRatioAndAssertInteger(relativeOffset.y,
                                                                                   proxy->dx),
                                                   tools::getRatioAndAssertInteger(relativeOffset.z,
                                                                                   proxy->dx),
                                                   target.size);
      }


      return proxy;
    }


    virtual void debugInfo(std::ostream &s) const {
      s << "Grid of side " << size << " address " << this << "; " << this->flags.size() << " cells marked";
    }

    virtual void getFlaggedCells(std::vector<size_t> &targetArray) const {
      targetArray.insert(targetArray.end(), flags.begin(), flags.end());
    }

    virtual void flagCells(const std::vector<size_t> &sourceArray) {
      // DEBUG:
      for(size_t x: sourceArray) {
        assert(x<size3);
      }
      std::vector<size_t> newFlags(flags.size() + sourceArray.size());
      auto end = std::set_union(flags.begin(), flags.end(), sourceArray.begin(), sourceArray.end(), newFlags.begin());
      newFlags.resize(end - newFlags.begin());
      flags = std::move(newFlags);
    }

    virtual void unflagAllCells() {
      flags.clear();
    }

    virtual bool hasFlaggedCells() const {
      return flags.size() > 0;
    }


    virtual bool pointsToGrid(const Grid *pOther) const {
      return this == pOther;
    }

    bool pointsToAnyGrid(std::vector<std::shared_ptr<Grid<T>>> grids) {
      for (auto g: grids) {
        if (pointsToGrid(g.get()))
          return true;
      }
      return false;
    }

    virtual bool containsPoint(const Coordinate<T> &coord) const {
      return Window<T>(simsize,offsetLower,offsetLower+boxsize).contains(coord);
    }

    virtual bool containsCellWithCoordinate(Coordinate<int> coord) const {
      return coord.x >= 0 && coord.y >= 0 && coord.z >= 0 &&
             coord.x < size && coord.y < size && coord.z < size;
    }

    virtual bool containsCell(size_t i) const {
      return i < size3;
    }

    virtual std::shared_ptr<Grid<T>> makeScaledMassVersion(T massRatio) {
      return std::make_shared<MassScaledGrid<T>>(this->shared_from_this(), massRatio);
    }


    ///////////////////////////////////////
    //  Index manipulation routines
    ///////////////////////////////////////

    size_t getIndexFromIndexAndStep(size_t index, const Coordinate<int> &step) const {
      auto coord = getCellCoordinate(index);
      coord += step;
      return this->getCellIndex(coord); // N.B. does wrapping inside getIndex
    }

    Coordinate<T> wrapPoint(Coordinate<T> pos) const {
      // TODO: this doesn't belong in the Grid class
      pos.x = fmod(pos.x, simsize);
      if (pos.x < 0) pos.x += simsize;
      pos.y = fmod(pos.y, simsize);
      if (pos.y < 0) pos.y += simsize;
      pos.z = fmod(pos.z, simsize);
      if (pos.z < 0) pos.z += simsize;
      return pos;
    }

    Coordinate<int> wrapCoordinate(Coordinate<int> index) const {
#ifdef SAFER_SLOWER
      index.x = index.x%simEquivalentSize;
      index.y = index.y%simEquivalentSize;
      index.z = index.z%simEquivalentSize;
#else
      if (index.x > (signed) simEquivalentSize - 1) index.x -= simEquivalentSize;
      if (index.y > (signed) simEquivalentSize - 1) index.y -= simEquivalentSize;
      if (index.z > (signed) simEquivalentSize - 1) index.z -= simEquivalentSize;
#endif
      if (index.x < 0) index.x += simEquivalentSize;
      if (index.y < 0) index.y += simEquivalentSize;
      if (index.z < 0) index.z += simEquivalentSize;
      return index;
    }


    size_t getCellIndex(Coordinate<int> coord) const {
      coord = wrapCoordinate(coord);
      return getCellIndexNoWrap(coord);
    }

    size_t getCellIndexNoWrap(size_t x, size_t y, size_t z) const {
      return (x * size + y) * size + z;
    }

    size_t getCellIndexNoWrap(int x, int y, int z) const {

#ifdef SAFER_SLOWER
      if(x<0 || x>=size || y<0 || y>=size || z<0 || z>=size)
          throw std::runtime_error("Grid index out of range in getIndexNoWrap");
#endif
      return size_t(x * size + y) * size + z;
    }

    size_t getCellIndexNoWrap(const Coordinate<int> &coordinate) const {
      return getCellIndexNoWrap(coordinate.x, coordinate.y, coordinate.z);
    }


    Coordinate<int> getCellCoordinate(size_t id) const {
      size_t x, y;

      if (id >= size3) {
        throw std::runtime_error("Index out of range");
      }

      // The following implementation is a little faster than using the
      // modulo operator.
      x = id / size2;
      id -= x * size2;
      y = id / size;
      id -= y * size;

      return Coordinate<int>(x, y, id);
    }

    Coordinate<int> getFourierCellCoordinate(size_t id) const {

      auto coord = getCellCoordinate(id);
      if (coord.x > (signed) size / 2) coord.x -= size;
      if (coord.y > (signed) size / 2) coord.y -= size;
      if (coord.z > (signed) size / 2) coord.z -= size;

      return coord;
    }

    T getFourierCellKSquared(size_t id) const {
      T res;
      auto coords = getFourierCellCoordinate(id);
      res = coords.x * coords.x + coords.y * coords.y + coords.z * coords.z;
      res *= kMinSquared;
      return res;
    }

    T getFourierCellAbsK(size_t id) const {
      return sqrt(getFourierCellKSquared(id));
    }

    T getFourierKmin() const {
      return kMin;
    }


    Coordinate<T> getCellCentroid(size_t id) const {
      Coordinate<int> coord = getCellCoordinate(id);
      return getCellCentroid(coord);
    }

    Coordinate<T> getCellCentroid(const Coordinate<int> & coord) const {
      Coordinate<T> result(coord);
      result *= dx;
      result += offsetLower;
      result += dx / 2;
      return result;
    }

    size_t getClosestId(Coordinate<T> coord) {
      auto coords = floor(wrapPoint(coord - offsetLower - dx / 2) / dx);
      return getCellIndexNoWrap(coords);
    }

    void appendIdsInCubeToVector(T x0c, T y0c, T z0c, T dxc, vector<size_t> &ids) {
      size_t offset = ids.size();
      int added_size = std::round(dxc/dx);
      added_size*=added_size*added_size;
      ids.resize(offset+added_size);
      insertCubeIdsIntoVector(x0c, y0c, z0c, dxc, ids.begin()+offset);
    }

    void insertCubeIdsIntoVector(T x0c, T y0c, T z0c, T dxc, vector<size_t>::iterator start) {
      // return all the grid IDs whose centres lie within the specified cube

      // TODO: optimization, set the storage size of ids here.

      std::tie(x0c, y0c, z0c) = wrapPoint(Coordinate<T>(x0c, y0c, z0c) - offsetLower);

      int xa = ((int) floor((x0c - dxc / 2 + dx / 2) / dx));
      int ya = ((int) floor((y0c - dxc / 2 + dx / 2) / dx));
      int za = ((int) floor((z0c - dxc / 2 + dx / 2) / dx));

      int xb = ((int) floor((x0c + dxc / 2 - dx / 2) / dx));
      int yb = ((int) floor((y0c + dxc / 2 - dx / 2) / dx));
      int zb = ((int) floor((z0c + dxc / 2 - dx / 2) / dx));


      iterateOverCube<int>(Coordinate<int>(xa, ya, za),
                           Coordinate<int>(xb, yb, zb) + 1,
                           [&start, this](const Coordinate<int> &cellCoord) {
                             (*start) = getCellIndex(cellCoord);
                             assert(*start<size3);
                             ++start;
                           });

    }

    int getEffectiveSimulationSize() const {
      return tools::getRatioAndAssertInteger(simsize, dx);
    }


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
        auto coord = source->getCellCoordinate(id);
        iterateOverCube<int>(
          coord * factor, coord * factor + factor,
          [&targetArray, &target](const Coordinate<int> &subCoord) {
            targetArray.push_back(target->getCellIndexNoWrap(subCoord));
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

      size_t factor = tools::getRatioAndAssertPositiveInteger(target->dx, source->dx);

      // it's not clear that the following parallelisation actually speeds much up
#pragma omp parallel for
      for (size_t i=0; i<sourceArray.size(); ++i) {
        size_t id = sourceArray[i];
        auto coord = source->getCellCoordinate(id);
        targetArray[i] = target->getCellIndexNoWrap(coord / factor);
      }

      // this sort seems to be the slowest step. In C++17 we can make it parallel... or is there a
      // better overall algorithm?
      std::sort(targetArray.begin(),targetArray.end());
      targetArray.erase(std::unique(targetArray.begin(),targetArray.end()),targetArray.end());
    }


  };

}

#endif

