#ifndef IC_VIRTUALGRID_HPP
#define IC_VIRTUALGRID_HPP

#include <cassert>
#include <set>
#include <type_traits>
#include <memory>
#include <vector>
#include <complex>

#include "src/simulation/coordinate.hpp"
#include "grid.hpp"
#include "src/simulation/window.hpp"

using std::complex;
using std::vector;
using std::make_shared;

namespace grids {
  template<typename T>
  class Grid;


  template<typename T>
  class VirtualGrid : public Grid<T> {
  protected:
    using typename Grid<T>::GridPtrType;
    using typename Grid<T>::ConstGridPtrType;
    GridPtrType pUnderlying;

  public:
    explicit VirtualGrid(GridPtrType pUnderlying) :
        Grid<T>(
            pUnderlying->periodicDomainSize, pUnderlying->size,
            pUnderlying->cellSize, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
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

    void debugInfo(std::ostream &s) const override {
      debugName(s);
      s << " of side " << this->size << " address " << this << " referencing ";
      pUnderlying->debugInfo(s);
    }

    bool pointsToGrid(const Grid<T> *pOther) const override {
      return this == pOther || pUnderlying.get() == pOther || pUnderlying->pointsToGrid(pOther);
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

    ConstGridPtrType getUnderlying() const {
      return pUnderlying;
    }


  };

  template<typename T>
  class SuperSampleGrid : public VirtualGrid<T> {
  private:
    int factor;
    int factor3;

  protected:
    using typename Grid<T>::GridPtrType;

  public:
    SuperSampleGrid(GridPtrType pUnderlying, int factor) :
        VirtualGrid<T>(pUnderlying,
                       pUnderlying->periodicDomainSize, pUnderlying->size * factor,
                       pUnderlying->cellSize / factor, pUnderlying->offsetLower.x,
                       pUnderlying->offsetLower.y,
                       pUnderlying->offsetLower.z,
                       pUnderlying->cellMassFrac / (factor * factor * factor),
                       factor),
        factor(factor) {

      factor3 = factor * factor * factor;
    }


    void debugName(std::ostream &s) const override {
      s << "SuperSampleGrid";
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray, upscaledArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      Grid<T>::upscaleCellFlagVector(underlyingArray, upscaledArray, this->pUnderlying.get(), this);
      targetArray.insert(targetArray.end(), upscaledArray.begin(), upscaledArray.end());
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::downscaleCellFlagVector(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->flagCells(targetArray);
    }


  };

  template<typename T>
  class ResolutionMatchingGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::GridPtrType;
    using typename Grid<T>::ConstGridPtrType;

    GridPtrType pUnderlyingHiRes, pUnderlyingLoResInterpolated;
    Coordinate<int> windowLowerCornerInclusive;
    Coordinate<int> windowUpperCornerExclusive;
    Coordinate<T> windowLowerCorner;
    Coordinate<T> windowUpperCorner;

  public:

    ResolutionMatchingGrid(GridPtrType pUnderlyingHiRes, GridPtrType pUnderlyingLoRes) :
        VirtualGrid<T>(pUnderlyingHiRes,
                       pUnderlyingLoRes->periodicDomainSize,
                       pUnderlyingLoRes->size *
                       tools::getRatioAndAssertInteger(pUnderlyingLoRes->cellSize, pUnderlyingHiRes->cellSize),
                       pUnderlyingHiRes->cellSize,
                       pUnderlyingLoRes->offsetLower.x,
                       pUnderlyingLoRes->offsetLower.y,
                       pUnderlyingLoRes->offsetLower.z,
                       pUnderlyingHiRes->cellMassFrac,
                       pUnderlyingHiRes->cellSofteningScale),
        pUnderlyingHiRes(pUnderlyingHiRes) {
      pUnderlyingLoResInterpolated = std::make_shared<SuperSampleGrid<T>>(pUnderlyingLoRes,
                                                                          this->size / pUnderlyingLoRes->size);
      this->pUnderlying = pUnderlyingLoResInterpolated;

      auto offsetLowerRelative = pUnderlyingHiRes->offsetLower - pUnderlyingLoRes->offsetLower;
      windowLowerCornerInclusive = round<int>(offsetLowerRelative / this->cellSize);
      windowLowerCorner = this->getCentroidFromCoordinate(windowLowerCornerInclusive) - this->cellSize / 2;

      auto offsetLowerRelativeCheck = Coordinate<T>(windowLowerCornerInclusive) * this->cellSize;
      assert(
          offsetLowerRelativeCheck.almostEqual(offsetLowerRelative)); // if this doesn't match, the grids don't line up

      windowUpperCornerExclusive = windowLowerCornerInclusive + pUnderlyingHiRes->size;
      windowUpperCorner = this->getCentroidFromCoordinate(windowUpperCornerExclusive) - this->cellSize / 2;

    }

    ConstGridPtrType getUnderlyingLoResInterpolated() const {
      return pUnderlyingLoResInterpolated;
    }

    ConstGridPtrType getUnderlyingHiRes() const {
      return pUnderlyingHiRes;
    }

    bool pointsToGrid(const Grid<T> *pOther) const override {
      return this == pOther || pUnderlyingLoResInterpolated->pointsToGrid(pOther) ||
             pUnderlyingHiRes->pointsToGrid(pOther);
    }

    void debugName(std::ostream &s) const override {
      s << "ResolutionMatchingGrid";
    }

    void debugInfo(std::ostream &s) const override {
      debugName(s);
      s << " of side " << this->size << " address " << this << " referencing (for genuine hi-res part) ";
      pUnderlyingHiRes->debugInfo(s);
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> interpolatedCellsArray;
      this->pUnderlyingHiRes->getFlaggedCells(targetArray);
      this->pUnderlyingLoResInterpolated->getFlaggedCells(interpolatedCellsArray);
      for (size_t i = 0; i < targetArray.size(); ++i) {
        auto coordinate =
            this->pUnderlyingHiRes->getCoordinateFromIndex(targetArray[i]) + this->windowLowerCornerInclusive;
        targetArray[i] = this->pUnderlyingLoResInterpolated->getIndexFromCoordinate(coordinate);
      }

      for (size_t i = 0; i < interpolatedCellsArray.size(); ++i) {
        size_t thisCell = interpolatedCellsArray[i];
        auto coordinate = pUnderlyingLoResInterpolated->getCoordinateFromIndex(thisCell);
        if (!coordinate.inWindow(windowLowerCornerInclusive, windowUpperCornerExclusive)) {
          targetArray.push_back(thisCell);
        }
      }

      tools::sortAndEraseDuplicate(targetArray);
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> interpolatedCellsArray;
      std::vector<size_t> hiresCellsArray;

      for (size_t i = 0; i < sourceArray.size(); ++i) {
        size_t thisCell = sourceArray[i];
        auto coordinate = pUnderlyingLoResInterpolated->getCoordinateFromIndex(thisCell);
        if (isInHiResWindow(coordinate)) {
          size_t hiresCell = getIndexInHiResWindow(coordinate);
          hiresCellsArray.push_back(hiresCell);
        } else {
          interpolatedCellsArray.push_back(thisCell);
        }
      }

      pUnderlyingLoResInterpolated->flagCells(interpolatedCellsArray);
      pUnderlyingHiRes->flagCells(hiresCellsArray);
    }

    bool isInHiResWindow(const Coordinate<int> &coordinate) const {
      return coordinate.inWindow(windowLowerCornerInclusive, windowUpperCornerExclusive);
    }

    bool isInHiResWindow(const Coordinate<T> &location) const {
      return location.inWindow(windowLowerCorner, windowUpperCorner);
    }

    size_t getIndexInHiResWindow(const Coordinate<int> &coordinate) const {
      return pUnderlyingHiRes->getIndexFromCoordinate(coordinate - windowLowerCornerInclusive);
    }


  };


  template<typename T>
  class SectionOfGrid : public VirtualGrid<T> {
  private:
    Coordinate<int> cellOffset;
    Coordinate<int> upperCell;
    Coordinate<T> posOffset;

  protected:
    using typename Grid<T>::GridPtrType;


  public:
    size_t mapIndexToUnderlying(size_t sec_id) const {
      auto coord = this->getCoordinateFromIndex(sec_id);
      coord += cellOffset;
      coord = this->wrapCoordinate(coord);
      if (!this->pUnderlying->containsCellWithCoordinate(coord))
        throw std::out_of_range("Out of range in SectionOfGrid::mapIndexToUnderlying");
      return this->pUnderlying->getIndexFromCoordinate(coord);
    }

    size_t mapIndexFromUnderlying(size_t underlying_id) const {
      Coordinate<int> coord = this->pUnderlying->getCoordinateFromIndex(underlying_id);
      coord -= cellOffset;
      coord = this->wrapCoordinate(coord);

      if (!this->containsCellWithCoordinate(coord))
        throw std::out_of_range("Out of range in SectionOfGrid::mapIndexFromUnderlying");

      return this->getIndexFromCoordinate(coord);
    }

    SectionOfGrid(GridPtrType pUnderlying, int deltax, int deltay, int deltaz, size_t size) :
        VirtualGrid<T>(pUnderlying,
                       pUnderlying->periodicDomainSize, size,
                       pUnderlying->cellSize,
                       pUnderlying->offsetLower.x + deltax * pUnderlying->cellSize,
                       pUnderlying->offsetLower.y + deltay * pUnderlying->cellSize,
                       pUnderlying->offsetLower.z + deltaz * pUnderlying->cellSize,
                       pUnderlying->cellMassFrac, pUnderlying->cellSofteningScale),
        cellOffset(deltax, deltay, deltaz),
        upperCell(cellOffset + pUnderlying->size),
        posOffset(deltax * this->cellSize, deltay * this->cellSize, deltaz * this->cellSize) {

    }


    bool containsPoint(const Coordinate<T> &coord) const override {
      return VirtualGrid<T>::containsPoint(coord) && this->pUnderlying->containsPoint(coord);
    }

    virtual bool containsCellWithCoordinate(const Coordinate<int> &coord) const {
      return VirtualGrid<T>::containsCellWithCoordinate(coord) &&
             this->pUnderlying->containsCellWithCoordinate(this->wrapCoordinate(coord + cellOffset));
    }

    virtual bool containsCell(size_t i) const {
      auto coord = this->getCoordinateFromIndex(i);
      return containsCellWithCoordinate(coord);
    }


    void debugName(std::ostream &s) const override {
      s << "SectionOfGrid";
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      for (size_t ptcl: underlyingArray) {
        try {
          targetArray.push_back(this->mapIndexFromUnderlying(ptcl));
          assert(targetArray.back() < this->size3);
        } catch (std::out_of_range &e) {
          continue;
        }
      }
     tools::sortAndEraseDuplicate(targetArray);
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> underlyingArray;
      for (size_t ptcl: sourceArray) {
        try {
          underlyingArray.push_back(this->mapIndexToUnderlying(ptcl));
        } catch (std::out_of_range &e) {
          continue;
        }
      }
      tools::sortAndEraseDuplicate(underlyingArray);
      this->pUnderlying->flagCells(underlyingArray);
    }


  };


  template<typename T>
  class SubSampleGrid : public VirtualGrid<T> {
  private:
    int factor;
    int factor3;

  public:
    SubSampleGrid(std::shared_ptr<Grid<T>> pUnderlying, int factor) :
        VirtualGrid<T>(pUnderlying,
                       pUnderlying->periodicDomainSize, pUnderlying->size / factor,
                       pUnderlying->cellSize * factor, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                       pUnderlying->offsetLower.z, pUnderlying->cellMassFrac * powf(factor, 3.0),
                       pUnderlying->cellSofteningScale),
        factor(factor) {
      //if(this->pUnderlying->size%factor!=0)
      //  throw std::runtime_error("SubSampleGrid - factor must be a divisor of the original grid size");

      factor3 = factor * factor * factor;
    }

   void debugName(std::ostream &s) const override {
      s << "SubSampleGrid";
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray, downscaledArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      Grid<T>::downscaleCellFlagVector(underlyingArray, downscaledArray, this->pUnderlying.get(), this);
      targetArray.insert(targetArray.end(), downscaledArray.begin(), downscaledArray.end());
      // err << "SubSample getFlaggedCells - underlying = " << underlyingArray.size() << " transformed = " <<targetArray.size() << endl;
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::upscaleCellFlagVector(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->flagCells(targetArray);
      // cerr << "SubSample flagCells - source = " << sourceArray.size() << " transformed = " <<targetArray.size() << endl;
    }


    int forEachSubcell(size_t id, std::function<void(size_t)> callback) const {
      auto coord0 = this->getCoordinateFromIndex(id);
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
            callback(this->pUnderlying->getIndexFromCoordinateNoWrap(xi, yi, zi));
          }

        }
      }
      return localFactor3;
    }


  };


  template<typename T>
  class MassScaledGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::GridPtrType;


  public:
    MassScaledGrid(GridPtrType pUnderlying, T massScale) :
        VirtualGrid<T>(pUnderlying,
                       pUnderlying->periodicDomainSize, pUnderlying->size,
                       pUnderlying->cellSize, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                       pUnderlying->offsetLower.z, massScale * pUnderlying->cellMassFrac,
                       pUnderlying->cellSofteningScale) {}

    void debugName(std::ostream &s) const override {
      s << "MassScaledGrid";
    }
  };

  //! Wrap a grid such that its center is a given point.
  /*! Does not change the offsetLower of the grid, i.e. do not move the grid wrto to its parent.
   * This class wraps coordinate and centroids but does not change the id of the cell, i.e.
   *  a cell (id, coord, centroid) in the underlying grid is mapped to (id, centered_coord, centered_centroid)
   *  in the virtual grid.
   */
  template<typename T>
  class CenteredGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::GridPtrType;

  private:
    const Coordinate<int> offset;
  public:
    CenteredGrid(GridPtrType pUnderlying, Coordinate<T> center) :
        VirtualGrid<T>(pUnderlying,
                       pUnderlying->periodicDomainSize, pUnderlying->size,
                       pUnderlying->cellSize,
                       pUnderlying->offsetLower.x,
                       pUnderlying->offsetLower.y,
                       pUnderlying->offsetLower.z,
                       pUnderlying->cellMassFrac,
                       pUnderlying->cellSofteningScale),
        offset((Coordinate<T>(0.5 * this->pUnderlying->thisGridSize) - center) / this->pUnderlying->cellSize ){}

    void debugName(std::ostream &s) const override {
      s << "CenteredGrid";
    }

    Coordinate<T> getPointOffset() const {
      return this->getPointOffsetfromCoordinateOffset(this->offset);
    }

  private:
    //! Transforming int offset in floating point offset
    Coordinate<T> getPointOffsetfromCoordinateOffset(Coordinate<int> offset) const {
      return Coordinate<T>(offset.x * this->pUnderlying->cellSize,
                           offset.y * this->pUnderlying->cellSize,
                           offset.z * this->pUnderlying->cellSize);
    }

    //! Get inverse centering transformation
    Coordinate<int> getInverseOffset() const {
      return Coordinate<int>(- this->offset.x, - this->offset.y, - this->offset.z);
    }

    Coordinate<T> getInversePointOffset() const {
      return this->getPointOffsetfromCoordinateOffset(this->getInverseOffset());
    }

  protected:

    size_t getIndexFromCoordinateNoWrap(size_t x, size_t y, size_t z) const override{
      return this->pUnderlying->getIndexFromIndexAndStep(
          this->pUnderlying->getIndexFromCoordinateNoWrap(x,y,z), this->getInverseOffset());
    }

    size_t getIndexFromCoordinateNoWrap(int x, int y, int z) const override {
      return this->pUnderlying->getIndexFromIndexAndStep(
          this->pUnderlying->getIndexFromCoordinateNoWrap(x,y,z), this->getInverseOffset());
    }

    size_t getIndexFromCoordinateNoWrap(const Coordinate<int> &coordinate) const override {
      return this->pUnderlying->getIndexFromIndexAndStep(
          this->pUnderlying->getIndexFromCoordinateNoWrap(coordinate), this->getInverseOffset());
    }

    size_t getIndexFromCoordinate(Coordinate<int> coord) const override{
      coord = this->pUnderlying->wrapCoordinate(coord);
      return this->getIndexFromCoordinateNoWrap(coord);
    }

    Coordinate<int> getCoordinateFromIndex(size_t id) const override {
      return this->pUnderlying->wrapCoordinate(this->pUnderlying->getCoordinateFromIndex(id) + this->offset);
    }

    Coordinate<T> getCentroidFromCoordinate(const Coordinate<int> &coord) const override {
      return this->pUnderlying->getCentroidFromCoordinate(
          this->pUnderlying->wrapCoordinate(this->offset + coord));
    }

    Coordinate<T> getCentroidFromIndex(size_t id) const override{
      return this->wrapPoint(this->pUnderlying->getCentroidFromIndex(id) + this->getPointOffset());
    }

    size_t getIndexFromPoint(Coordinate<T> point) const override{
      return this->pUnderlying->getIndexFromPoint(point + this->getInversePointOffset());
    }

    Coordinate<int> getCoordinateFromPoint(Coordinate<T> point) const override{
      return this->pUnderlying->getCoordinateFromPoint(point);
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      this->pUnderlying->getFlaggedCells(targetArray);
    }
  };

  //! Offsets the corner of a grid in the cordinate of the box
  /*! Modifies only the centroid coordinates, e.g. a cell (id, coord, centroid) in the underlying grid
   * will be mapped to (id, coord, centroid + offset) by te virtual grid.
   */
  template<typename T>
  class OffsetGrid : public VirtualGrid<T> {

  protected:
    using typename Grid<T>::GridPtrType;


  public:
    OffsetGrid(GridPtrType pUnderlying, T dx, T dy, T dz) :
        VirtualGrid<T>(pUnderlying, pUnderlying->periodicDomainSize, pUnderlying->size,
                       pUnderlying->cellSize,
                       pUnderlying->offsetLower.x + dx,
                       pUnderlying->offsetLower.y + dy,
                       pUnderlying->offsetLower.z + dz,
                       pUnderlying->cellMassFrac,
                       pUnderlying->cellSofteningScale) {}


    void debugName(std::ostream &s) const override {
      s << "OffsetGrid";
    }
  };




}


#endif
