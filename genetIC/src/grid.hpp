#ifndef __GRID_HPP
#define __GRID_HPP

#include <cassert>
#include <set>
#include <type_traits>
#include <memory>
#include <vector>
#include <complex>
#include "src/numerics/fourier.hpp"
#include "filter.hpp"
#include "coordinate.hpp"
#include "progress/progress.hpp"
#include "utils.hpp"

using std::complex;
using std::vector;
using std::make_shared;

namespace fields {
  template<typename T, typename S>
  class Field;

  template<typename T, typename S>
  class InterpolatedField;

  template<typename T, typename S>
  class SubsampledField;

  template<typename T, typename S>
  class TranslatedField;
}

namespace grids {
  template<typename T>
  class Grid;

  template<typename T>
  class SuperSampleGrid;

  template<typename T>
  class SubSampleGrid;

  template<typename T>
  class OffsetGrid;

  template<typename T>
  class SectionOfGrid;

  template<typename T>
  class MassScaledGrid;


  template<typename T>
  class Grid : public std::enable_shared_from_this<Grid<T>> {
  public:

    using TField = ::fields::Field<std::complex<T>, T>;
    using TRealField = ::fields::Field<T, T>;
    using PtrTField = std::shared_ptr<std::vector<std::complex<T>>>;
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
    const T cellMassFrac; ///< the fraction of mass of the full simulation in a single cell of this grid
    const T cellSofteningScale; ///< normally 1.0; scales softening relative to dx

    Grid(T simsize, size_t n, T dx = 1.0, T x0 = 0.0, T y0 = 0.0, T z0 = 0.0,
         T massFrac = 0.0, T softScale = 1.0) :
      simsize(simsize), boxsize(dx * n),
      dx(dx), offsetLower(x0, y0, z0),
      size(n), size2(n * n), size3(n * n * n),
      cellMassFrac(massFrac == 0.0 ? pow(dx / simsize, 3.0) : massFrac),
      cellSofteningScale(softScale) {
      // cerr << "Grid ctor " << this <<  endl;
      setKmin();
    }


    Grid(size_t n) : simsize(0), boxsize(n),
                     dx(1.0), offsetLower(0, 0, 0),
                     size(n), size2(n * n), size3(n * n * n), cellMassFrac(0.0),
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

    GridPtrType makeProxyGridToMatch(const Grid<T> &target) const {
      GridPtrType proxy = std::const_pointer_cast<Grid<T>>(this->shared_from_this());
      if (target.dx > dx) {
        size_t ratio = getRatioAndAssertPositiveInteger(target.dx, dx);
        proxy = std::make_shared<SubSampleGrid<T>>(proxy, ratio);
      } else if (target.dx < dx) {
        size_t ratio = getRatioAndAssertPositiveInteger(dx, target.dx);
        proxy = std::make_shared<SuperSampleGrid<T>>(proxy, ratio);
      }

      if (target.offsetLower != offsetLower || target.size != proxy->size) {
        proxy = std::make_shared<SectionOfGrid<T>>(proxy,
                                                   getRatioAndAssertInteger(target.offsetLower.x - offsetLower.x,
                                                                            proxy->dx),
                                                   getRatioAndAssertInteger(target.offsetLower.y - offsetLower.y,
                                                                            proxy->dx),
                                                   getRatioAndAssertInteger(target.offsetLower.z - offsetLower.z,
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
      flags.insert(flags.end(), sourceArray.begin(), sourceArray.end());
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

    virtual bool containsCell(const Coordinate<T> &coord) const {
      return coord.x >= 0 && coord.y >= 0 && coord.z >= 0 &&
             coord.x < size && coord.y < size && coord.z < size;
    }

    virtual bool containsCell(size_t i) const {
      return i < size3;
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) const {
      return field[i];
    }

    virtual T getFieldAt(size_t i, const TRealField &field) const {
      return field[i];
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

    void simWrap(Coordinate<T> &pos) const {
      // TODO: this doesn't belong in the Grid class
      pos.x = fmod(pos.x, simsize);
      if (pos.x < 0) pos.x += simsize;
      pos.y = fmod(pos.y, simsize);
      if (pos.y < 0) pos.y += simsize;
      pos.z = fmod(pos.z, simsize);
      if (pos.z < 0) pos.z += simsize;
    }

    void cellWrap(Coordinate<int> &index) const {
#ifdef SAFER_SLOWER
      index.x = index.x%size;
      index.y = index.y%size;
      index.z = index.z%size;
#else
      if (index.x > (signed) size - 1) index.x -= size;
      if (index.y > (signed) size - 1) index.y -= size;
      if (index.z > (signed) size - 1) index.z -= size;
#endif
      if (index.x < 0) index.x += size;
      if (index.y < 0) index.y += size;
      if (index.z < 0) index.z += size;
    }


    size_t getCellIndex(Coordinate<int> coord) const {
      cellWrap(coord);
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
      int x, y, z;

      if (id >= size3) throw std::runtime_error("Index out of range");

      // The following implementation is a little faster than using the
      // modulo operator.
      x = int(id / size2);
      id -= size_t(x) * size2;
      y = int(id / size);
      id -= size_t(y) * size;
      z = int(id);

      return Coordinate<int>(x, y, z);
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


    Coordinate<T> getCellCentroid(size_t id) const {
      Coordinate<T> coord = getCellCoordinate(id);
      coord *= dx;
      coord += offsetLower;
      coord += dx / 2;
      return coord;
    }

    size_t getClosestIdNoWrap(Coordinate<T> coord) {
      auto coords = floor((coord - offsetLower - dx / 2) / dx);
      return getCellIndexNoWrap(coords);
    }

    void appendIdsInCubeToVector(T x0c, T y0c, T z0c, T dxc, vector<size_t> &ids) {
      // return all the grid IDs whose centres lie within the specified cube

      // TODO: optimization, set the storage size of ids here.

      int xa = ((int) floor((x0c - offsetLower.x - dxc / 2 + dx / 2) / dx));
      int ya = ((int) floor((y0c - offsetLower.y - dxc / 2 + dx / 2) / dx));
      int za = ((int) floor((z0c - offsetLower.z - dxc / 2 + dx / 2) / dx));

      int xb = ((int) floor((x0c - offsetLower.x + dxc / 2 - dx / 2) / dx));
      int yb = ((int) floor((y0c - offsetLower.y + dxc / 2 - dx / 2) / dx));
      int zb = ((int) floor((z0c - offsetLower.z + dxc / 2 - dx / 2) / dx));

      iterateOverCube<int>(Coordinate<int>(xa, ya, za),
                           Coordinate<int>(xb, yb, zb) + 1,
                           [&ids, this](const Coordinate<int> &cellCoord) {
                             ids.emplace_back(getCellIndex(cellCoord));
                           });

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


      std::set<size_t> targetSet;

      assert(source->size >= target->size);
      assert(source->offsetLower == target->offsetLower);

      assert((source->size) % (target->size) == 0);
      size_t factor = source->size / target->size;

      for (auto id: sourceArray) {
        auto coord = source->getCellCoordinate(id);
        targetSet.insert(target->getCellIndexNoWrap(coord / factor));
      }
      targetArray.clear();
      targetArray.insert(targetArray.end(), targetSet.begin(), targetSet.end());
    }


  };


  template<typename T>
  class VirtualGrid : public Grid<T> {
  protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;
    GridPtrType pUnderlying;

  public:
    VirtualGrid(GridPtrType pUnderlying) :
      Grid<T>(
        pUnderlying->simsize, pUnderlying->size,
        pUnderlying->dx, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
        pUnderlying->offsetLower.z, pUnderlying->cellMassFrac, pUnderlying->cellSofteningScale),
      pUnderlying(pUnderlying) {

    }

    VirtualGrid(GridPtrType pUnderlying, T simsize, T gridsize,
                T dx, T x0, T y0, T z0, T massfrac = 0.0, T softScale = 1.0) :
      Grid<T>(simsize, gridsize, dx, x0, y0, z0, massfrac, softScale),
      pUnderlying(pUnderlying) {

    }

    virtual void debugName(std::ostream &s) const {
      s << "VirtualGrid";
    }

    virtual void debugInfo(std::ostream &s) const override {
      debugName(s);
      s << " of side " << this->size << " address " << this << " referencing ";
      pUnderlying->debugInfo(s);
    }

    bool pointsToGrid(const Grid<T> *pOther) const override {
      return pUnderlying.get() == pOther || pUnderlying->pointsToGrid(pOther);
    }

    virtual T getFieldAt(size_t i, const TRealField &field) const override {
      return pUnderlying->getFieldAt(i, field);
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) const override {
      return pUnderlying->getFieldAt(i, field);
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      pUnderlying->getFlaggedCells(targetArray);
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      pUnderlying->flagCells(sourceArray);
    }

    void unflagAllCells() override {
      pUnderlying->unflagAllCells();
    }

    bool hasFlaggedCells() const override {
      return pUnderlying->hasFlaggedCells();
    }


  };

  template<typename T>
  class SuperSampleGrid : public VirtualGrid<T> {
  private:
    int factor;
    int factor3;

  protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;

  public:
    SuperSampleGrid(GridPtrType pUnderlying, int factor) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->simsize, pUnderlying->size * factor,
                     pUnderlying->dx / factor, pUnderlying->offsetLower.x,
                     pUnderlying->offsetLower.y,
                     pUnderlying->offsetLower.z,
                     pUnderlying->cellMassFrac / (factor * factor * factor),
                     factor),
      factor(factor) {

      factor3 = factor * factor * factor;
    }


    virtual void debugName(std::ostream &s) const override {
      s << "SuperSampleGrid";
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      Grid<T>::upscaleCellFlagVector(underlyingArray, targetArray, this->pUnderlying.get(), this);
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::downscaleCellFlagVector(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->flagCells(targetArray);
    }

    virtual T getFieldAt(size_t i, const TRealField &field) const override {
      auto centroid = this->getCellCentroid(i);
      return field.evaluateInterpolated(centroid);
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) const override {
      auto centroid = this->getCellCentroid(i);
      return field.evaluateInterpolated(centroid);
    }

  };

  template<typename T>
  class OffsetGrid : public VirtualGrid<T> {

  protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;


  public:
    OffsetGrid(GridPtrType pUnderlying, T dx, T dy, T dz) :
      VirtualGrid<T>(pUnderlying, pUnderlying->simsize, pUnderlying->size,
                     pUnderlying->dx,
                     pUnderlying->offsetLower.x + dx,
                     pUnderlying->offsetLower.y + dy,
                     pUnderlying->offsetLower.z + dz,
                     pUnderlying->cellMassFrac,
                     pUnderlying->cellSofteningScale) {

    }


    virtual void debugName(std::ostream &s) const override {
      s << "OffsetGrid";
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      throw (std::runtime_error("getFlaggedCells is not implemented for OffsetGrid"));
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      throw (std::runtime_error("flagCells is not implemented for OffsetGrid"));
    }


  };


  template<typename T>
  class SectionOfGrid : public VirtualGrid<T> {
  private:
    Coordinate<int> cellOffset;
    Coordinate<int> upperCell;
    Coordinate<T> posOffset;

  protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;


    size_t mapIndexToUnderlying(size_t sec_id) const {
      auto coord = this->getCellCoordinate(sec_id);
      coord += cellOffset;
      if (!this->pUnderlying->containsCell(coord))
        throw std::out_of_range("Out of range in SectionOfGrid::mapIndexToUnderlying");
      return this->pUnderlying->getCellIndex(coord);
    }

    size_t mapIndexFromUnderlying(size_t underlying_id) const {
      auto coord = this->pUnderlying->getCellCoordinate(underlying_id);
      coord -= cellOffset;

      if (!this->containsCell(coord))
        throw std::out_of_range("Out of range in SectionOfGrid::mapIndexFromUnderlying");

      return this->getCellIndex(coord);
    }

  public:
    SectionOfGrid(GridPtrType pUnderlying, int deltax, int deltay, int deltaz, size_t size) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->simsize, size,
                     pUnderlying->dx,
                     pUnderlying->offsetLower.x + deltax * pUnderlying->dx,
                     pUnderlying->offsetLower.y + deltay * pUnderlying->dx,
                     pUnderlying->offsetLower.z + deltaz * pUnderlying->dx,
                     pUnderlying->cellMassFrac, pUnderlying->cellSofteningScale),
      cellOffset(deltax, deltay, deltaz),
      upperCell(cellOffset + pUnderlying->size),
      posOffset(deltax * this->dx, deltay * this->dx, deltaz * this->dx) {

    }


    virtual bool containsCell(size_t i) const override {
      return containsCell(this->getCellCoordinate(i));
    }

    virtual bool containsCell(const Coordinate<T> &coord) const override {
      auto translated_coord = coord + cellOffset;
      return this->pUnderlying->containsCell(translated_coord);
    }


    virtual void debugName(std::ostream &s) const override {
      s << "SectionOfGrid";
    }

    virtual T getFieldAt(size_t i, const TRealField &field) const override {
      return this->pUnderlying->getFieldAt(mapIndexToUnderlying(i), field);
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) const override {
      return this->pUnderlying->getFieldAt(mapIndexToUnderlying(i), field);
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      for (size_t ptcl: underlyingArray) {
        try {
          targetArray.push_back(this->mapIndexFromUnderlying(ptcl));
        } catch (std::out_of_range &e) {
          continue;
        }
      }
      std::sort(targetArray.begin(), targetArray.end());
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      throw (std::runtime_error("flagCells is not implemented for OffsetGrid"));
    }


  };


  template<typename T>
  class SubSampleGrid : public VirtualGrid<T> {
  private:
    int factor;
    int factor3;

  protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;

  public:
    SubSampleGrid(std::shared_ptr<Grid<T>> pUnderlying, int factor) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->simsize, pUnderlying->size / factor,
                     pUnderlying->dx * factor, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                     pUnderlying->offsetLower.z, pUnderlying->cellMassFrac * powf(factor, 3.0),
                     pUnderlying->cellSofteningScale),
      factor(factor) {
      //if(this->pUnderlying->size%factor!=0)
      //  throw std::runtime_error("SubSampleGrid - factor must be a divisor of the original grid size");

      factor3 = factor * factor * factor;
    }

    virtual void debugName(std::ostream &s) const override {
      s << "SubSampleGrid";
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      Grid<T>::downscaleCellFlagVector(underlyingArray, targetArray, this->pUnderlying.get(), this);
      // err << "SubSample getFlaggedCells - underlying = " << underlyingArray.size() << " transformed = " <<targetArray.size() << endl;
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::upscaleCellFlagVector(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->flagCells(targetArray);
      // cerr << "SubSample flagCells - source = " << sourceArray.size() << " transformed = " <<targetArray.size() << endl;
    }


    int forEachSubcell(size_t id, std::function<void(size_t)> callback) const {
      auto coord0 = this->getCellCoordinate(id);
      coord0 *= factor;
      auto coord1 = coord0 + factor;
      // In case the edge of the last cell in the fine grid (this->pUnderlying)
      // is not aligned with the edge of any cell in the coarse grid (this),
      // we need to be able to do an average over fewer than factor^3 cells.
      //
      // This is a situation we might, in retrospect, wish to avoid. However,
      // since early tests and simulations were conducted without obeying
      // this 'end-alignment' condition, we need to support it.
      const int underlyingSize = int(this->pUnderlying->size);
      if (coord1.x > underlyingSize) coord1.x = underlyingSize;
      if (coord1.y > underlyingSize) coord1.y = underlyingSize;
      if (coord1.z > underlyingSize) coord1.z = underlyingSize;

      int localFactor3 = (coord1.x - coord0.x) * (coord1.y - coord0.y) * (coord1.z - coord0.z);

      for (auto xi = coord0.x; xi < coord1.x; ++xi) {
        for (auto yi = coord0.y; yi < coord1.y; ++yi) {
          for (auto zi = coord0.z; zi < coord1.z; ++zi) {
            callback(this->pUnderlying->getCellIndexNoWrap(xi, yi, zi));
          }

        }
      }
      return localFactor3;
    }


    virtual complex<T> getFieldAt(size_t i, const TField &field) const override {
      complex<T> returnVal(0);
      int localFactor3 = forEachSubcell(i, [this, &returnVal, &field](size_t local_id) {
        returnVal += field[local_id];
      });
      return returnVal / T(localFactor3);
    }

    virtual T getFieldAt(size_t i, const TRealField &field) const override {
      T returnVal(0);
      int localFactor3 = forEachSubcell(i, [this, &returnVal, &field](size_t local_id) {
        returnVal += field[local_id];
      });
      return returnVal / localFactor3;
    }


  };


  template<typename T>
  class MassScaledGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;


  public:
    MassScaledGrid(GridPtrType pUnderlying, T massScale) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->simsize, pUnderlying->size,
                     pUnderlying->dx, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                     pUnderlying->offsetLower.z, massScale * pUnderlying->cellMassFrac,
                     pUnderlying->cellSofteningScale) {}

    virtual void debugName(std::ostream &s) const override {
      s << "MassScaledGrid";
    }


  };
}

#endif

