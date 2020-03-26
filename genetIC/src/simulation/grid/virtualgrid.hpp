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

  /*! \brief Class that defines a 'virtual' grid that sits on top of another one, and may have different properties.
    For example - virtual grids may describe the same grid with a higher resolution (super sampled grid), or with
    a lower resolution (sub-sampled grid), or simply with the centre of the grid offset by come vector.

  */
  template<typename T>
  class VirtualGrid : public Grid<T> {
  protected:
    using typename Grid<T>::GridPtrType;
    using typename Grid<T>::ConstGridPtrType;
    GridPtrType pUnderlying; //!< Pointer to the 'real' grid underlying this virtual grid.

  public:
    //! Construct a virtual grid with identical properties to the underlying grid.
    explicit VirtualGrid(GridPtrType pUnderlying) :
      Grid<T>(
        pUnderlying->periodicDomainSize, pUnderlying->size,
        pUnderlying->cellSize, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
        pUnderlying->offsetLower.z, pUnderlying->cellMassFrac, pUnderlying->cellSofteningScale),
      pUnderlying(pUnderlying) {

    }

    /*! \brief Construct a virtual grid with different properties to the underlying grid.

        \param pUnderlying - pointer to the underlying grid.
        \param simsize - size of the simulation in comoving units.
        \param gridsize - size of this grid.
        \param dx - size of one cell in comoving units
        \param x0 - x co-ordinate of the offset of the lower front left hand corner.
        \param y0 - y co-ordinate of the offset of the lower front left hand corner.
        \param z0 - z co-ordinate of the offset of the lower front left hand corner.
        \param massfrac - fraction of the total simulation mass in each cell (default value 0 instructs grid constructor to compute this).
        \param softscale - cell softening scale used by the grid.
    */
    VirtualGrid(GridPtrType pUnderlying, T simsize, T gridsize,
                T dx, T x0, T y0, T z0, T massfrac = 0.0, T softScale = 1.0) :
      Grid<T>(simsize, gridsize, dx, x0, y0, z0, massfrac, softScale),
      pUnderlying(pUnderlying) {

    }

    //! If the virtual grid has a different cell-size as the underlying grid, returns true, otherwise, return whether the underlying grid was up/down-sampled.
    virtual bool isUpsampledOrDownsampled() override {
      if (this->pUnderlying->cellSize != this->cellSize) {
        return true;
      } else {
        return this->pUnderlying->isUpsampledOrDownsampled();
      }
    }

    //! Output the name of this grid type. For debugging purposes only.
    virtual void debugName(std::ostream &s) const {
      s << "VirtualGrid";
    }

    //! Output debug information to the specified stream, along with some additional information including the underlying grid's debug information.
    void debugInfo(std::ostream &s) const override {
      debugName(s);
      s << " of side " << this->size << " address " << this << " referencing ";
      pUnderlying->debugInfo(s);
    }

    //! Returns true if this is a virtual version of the specified grid (possibly via some tower of other virtual grids).
    bool isProxyFor(const Grid<T> *pOther) const override {
      return this == pOther || pUnderlying.get() == pOther || pUnderlying->isProxyFor(pOther);
    }

    virtual //! Returns the list of flagged cells in the underlying grid.
    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      pUnderlying->getFlaggedCells(targetArray);
    }

    virtual //! Flags the specified cells in the underlying grid.
    void flagCells(const std::vector<size_t> &sourceArray) override {
      pUnderlying->flagCells(sourceArray);
    }

    virtual //! Unflag all flagged cells in the underlying grid.
    void unflagAllCells() override {
      pUnderlying->unflagAllCells();
    }

    virtual //! Returns the number of flagged cells in the underlying grid.
    size_t numFlaggedCells() const override {
      std::vector<size_t> t;
      getFlaggedCells(t);
      return t.size();
    }

    //! Returns a constant pointer to the underlying grid.
    ConstGridPtrType getUnderlying() const {
      return pUnderlying;
    }


  };

  /*! \brief Virtual grid with a higher resolution than the underlying grid.
    Specifically, this class can be used to flag/unflag cells as if we were dealing with a higher resolution grid, and
    the resulting flags will be downscaled/upscaled respectively to the underlying grid.

    Field data is obtained by interpolating the underlying low resolution grid.
  */
  template<typename T>
  class SuperSampleGrid : public VirtualGrid<T> {
  private:
    int factor; //!< Cells of the SuperSampleGrid are factor times smaller than those of the underlying grid.
    int factor3; //!< SuperSampleGrid has factor3 times as many cells as the underlying grid.

  protected:
    using typename Grid<T>::GridPtrType;

  public:
    //! \brief Constructor - specify the underlying grid, and how many times higher resolution the SuperSampleGrid should have
    /*! \param pUnderlying - pointer to the underlying grid.
        \param factor - how many times higher resolution than the underlying grid that the SuperSampleGrid should have.
    */
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


    //! Outputs the grid debug name. For debugging purposes only.
    void debugName(std::ostream &s) const override {
      s << "SuperSampleGrid";
    }

    //!\brief Returns the cell flags stored on the underlying grid, upscaling to give a vector of the corresponding flags as if they lived on the SuperSampleGrid
    /*! \param targetArray - vector to store flagged cells interpreted as if in the superSampleGrid
    */
    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray, upscaledArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      Grid<T>::upscaleCellFlagVector(underlyingArray, upscaledArray, this->pUnderlying.get(), this);
      targetArray.insert(targetArray.end(), upscaledArray.begin(), upscaledArray.end());
    }

    //! \brief Flag the specified cells, referred to as if they lived in the SuperSampleGrid.
    /*! \param sourceArray - vector of cells to flag, as if they were in the SuperSampleGrid
    */
    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::downscaleCellFlagVector(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->flagCells(targetArray);
    }

    GridPtrType makeSubsampled(size_t ratio) const override {
      // Special case: subsampling a supersampled grid can needlessly destroy accuracy, so avoid it!
      if (ratio > factor) {
        // relative to the underlying grid we will be subsampled
        size_t ratioToUnderlying = tools::getRatioAndAssertInteger(static_cast<float>(ratio),
                                                                   static_cast<float>(factor));
        return this->pUnderlying->makeSubsampled(ratioToUnderlying);
      } else if (ratio == factor) {
        // gives back the underlying grid
        return std::const_pointer_cast<Grid<T>>(this->pUnderlying);
      } else {
        // relative to the underlying grid we will still be supersampled
        size_t ratioToUnderlying = tools::getRatioAndAssertInteger(static_cast<float>(factor),
                                                                   static_cast<float>(ratio));
        return this->pUnderlying->makeSupersampled(ratioToUnderlying);
      }
    }

    GridPtrType makeSupersampled(size_t ratio) const override {
      return this->pUnderlying->makeSupersampled(ratio * factor);
    }

    bool containsCellWithCoordinate(const Coordinate<int> &coord) const override {
      return this->pUnderlying->containsCellWithCoordinate(coord/this->factor);
    }

    bool containsCell(size_t i) const override {
      auto coord = this->getCoordinateFromIndex(i);
      return this->containsCellWithCoordinate(coord);
    }

  };

  /*! \brief Virtual grid with two underlying grids - one high resolution, and one low resolution
    The low resolution grid is super-sampled up to match the resolution of the high resolution grid, which is then used as the underlying grid.
    The high resolution grid is treated as a high-resolution window that lies somewhere inside the low resolution grid, and we can use the
    resolution matching grid to reference points in the grid as if they were all at the higher resolution, even though data is only stored at
    high resolution within the high-res window.
  */
  template<typename T>
  class ResolutionMatchingGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::GridPtrType;
    using typename Grid<T>::ConstGridPtrType;

    GridPtrType pUnderlyingHiRes; //!< Pointer to the underlying high resolution grid.
    GridPtrType pUnderlyingLoRes; //!< Pointer to the underling low resolution (but bigger) grid
    GridPtrType pUnderlyingLoResInterpolated; //!< Pointers to the low resolution grid, interpolated up to the same resolution as the high resolution grid.
    Coordinate<int> windowLowerCornerInclusive; //!< Lower front left corner of the high resolution window region
    Window<int> windowInCoordinates; //!< Window for the high resolution region (integer version)
    Window<T> windowInSpace; //!< Window in co-moving co-ordinates.



  public:

    ResolutionMatchingGrid(GridPtrType pUnderlyingHiRes, GridPtrType pUnderlyingLoRes) :
      VirtualGrid<T>(
        pUnderlyingHiRes, //!< Underlying high-res grid is initially the underlying grid, but ultimately will be an super-sampled version of the low res grid.
        pUnderlyingLoRes->periodicDomainSize, //!< Uses the low resolution grid to set the simulation size
        pUnderlyingLoRes->size *
        tools::getRatioAndAssertInteger(pUnderlyingLoRes->cellSize,
                                        pUnderlyingHiRes->cellSize), //!< Number of cells that would be on one side if the whole simulation were are the resolution of the high-res grid
        pUnderlyingHiRes->cellSize,//!< Size of a single cell matches that of the high res grid
        pUnderlyingLoRes->offsetLower.x,//!< x co-ord of lower front left hand corner matches that of low-res grid.
        pUnderlyingLoRes->offsetLower.y,//!< y co-ord of lower front left hand corner matches that of low-res grid.
        pUnderlyingLoRes->offsetLower.z,//!< z co-ord of lower front left hand corner matches that of low-res grid.
        pUnderlyingHiRes->cellMassFrac, //!< Mass in each cell matches that of high res grid
        pUnderlyingHiRes->cellSofteningScale), //!< Cell softening scale matches that of high res grid
      pUnderlyingHiRes(pUnderlyingHiRes), pUnderlyingLoRes(pUnderlyingLoRes) {

      Coordinate<int> windowUpperCornerExclusive;
      Coordinate<T> windowLowerCornerInSpace;
      Coordinate<T> windowUpperCornerInSpace;

      // Create a super-sampled virtual grid from the low res grid, and set it to be the underlying grid:
      pUnderlyingLoResInterpolated =
        pUnderlyingLoRes->makeSupersampled(this->size / pUnderlyingLoRes->size)->withIndependentFlags();
      this->pUnderlying = pUnderlyingLoResInterpolated;

      // Set the location of the high res windowed region to the lower front left corner of the high res grid, relative to the low-res grid:
      auto offsetLowerRelative = pUnderlyingHiRes->offsetLower - pUnderlyingLoRes->offsetLower;
      windowLowerCornerInclusive = round<int>(offsetLowerRelative / this->cellSize);
      windowLowerCornerInSpace = this->getCentroidFromCoordinate(windowLowerCornerInclusive) - this->cellSize / 2;

      // Check that the grids line up correctly:
      auto offsetLowerRelativeCheck = Coordinate<T>(windowLowerCornerInclusive) * this->cellSize;
      assert(
        offsetLowerRelativeCheck.almostEqual(offsetLowerRelative)); // if this doesn't match, the grids don't line up

      // Get the location of the upper-back right corner of the high res window:
      windowUpperCornerExclusive = this->wrapCoordinate(windowLowerCornerInclusive + pUnderlyingHiRes->size);
      windowUpperCornerInSpace = this->getCentroidFromCoordinate(windowUpperCornerExclusive) - this->cellSize / 2;

      // Create the high resolution window in co-ordinate and co-moving space:
      windowInCoordinates = Window<int>(this->simEquivalentSize, windowLowerCornerInclusive,
                                        windowUpperCornerExclusive);
      windowInSpace = Window<T>(this->periodicDomainSize, windowLowerCornerInSpace, windowUpperCornerInSpace);


    }

    //! Returns a constant pointer to the super-sampled underlying low resolution grid
    ConstGridPtrType getUnderlyingLoResInterpolated() const {
      return pUnderlyingLoResInterpolated;
    }

    //! Returns a constant pointer to the high resolution grid.
    ConstGridPtrType getUnderlyingHiRes() const {
      return pUnderlyingHiRes;
    }

    //! Returns true if the high resolution grid is a proxy for the other grid
    bool isProxyFor(const Grid<T> *pOther) const override {
      // Note here that we count this grid as a proxy only for the high resolution underlying grid
      // This is crucial for the correct behaviour with subsampling; if it claims to be a proxy for the
      // low resolution grid, it might get further downsampled.
      return this == pOther || pUnderlyingHiRes->isProxyFor(pOther);
    }

    //! Debug information (for debugging purposes only)
    void debugName(std::ostream &s) const override {
      s << "ResolutionMatchingGrid";
    }

    //! Output more debug information, including information about the high resolution grid.
    void debugInfo(std::ostream &s) const override {
      debugName(s);
      s << " of side " << this->size << " address " << this << " referencing (for genuine hi-res part) ";
      pUnderlyingHiRes->debugInfo(s);
      s << "; (for interpolated part) ";
      pUnderlyingLoResInterpolated->debugInfo(s);
    }

    //! Stores the flagged cells in the targetArray vector as if they were defined on the super-sampled low resolution grid.
    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> interpolatedCellsArray;
      // Store all the flagged cells from the high resolution region:
      this->pUnderlyingHiRes->getFlaggedCells(targetArray);
      this->pUnderlyingLoResInterpolated->getFlaggedCells(interpolatedCellsArray);
      // Map the co-ordinates of the flagged cells from the high resolution region into co-ordinates in the super-sampled low resolution grid
      // and store the relevant indices:
      for (size_t i = 0; i < targetArray.size(); ++i) {
        auto coordinate =
          this->pUnderlyingHiRes->getCoordinateFromIndex(targetArray[i]) + this->windowLowerCornerInclusive;
        // The coordinate will be wrapped (if necessary) within getIndexFromCoordinate
        targetArray[i] = this->pUnderlyingLoResInterpolated->getIndexFromCoordinate(coordinate);
      }

      // Get any cells flagged that were not in the high resolution region, and store them:
      for (size_t i = 0; i < interpolatedCellsArray.size(); ++i) {
        size_t thisCell = interpolatedCellsArray[i];
        auto coordinate = pUnderlyingLoResInterpolated->getCoordinateFromIndex(thisCell);
        if (!isInHiResWindow(coordinate)) {
          targetArray.push_back(thisCell);
        }
      }

      tools::sortAndEraseDuplicate(targetArray);
    }

    //! Flags all the cells specified in sourceArray as if they were cells in the super-sampled low resolution grid, converting them to flags on the two underlying grids
    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> interpolatedCellsArray;
      std::vector<size_t> hiresCellsArray;

      // Flag cells in the appropriate grid according to whether or not they lie in the high resolution window:
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

    //! Unflag all cells on the underlying grids
    void unflagAllCells() override {
      pUnderlyingLoResInterpolated->unflagAllCells();
      pUnderlyingHiRes->unflagAllCells();
    }

    //! Return true if the specified co-ordinate lies in the high resolution window
    bool isInHiResWindow(const Coordinate<int> &coordinate) const {
      return windowInCoordinates.contains(coordinate);
      //coordinate.inWindow(windowLowerCornerInclusive, windowUpperCornerExclusive);
    }

    //! Return true if the specified co-moving position vector points to somewhere in the high resolution window
    bool isInHiResWindow(const Coordinate<T> &location) const {
      return windowInSpace.contains(location);
      //location.inWindow(windowLowerCorner, windowUpperCorner);
    }

    //! Returns the index of the specified co-ordinate relative to the high resolution grid
    size_t getIndexInHiResWindow(const Coordinate<int> &coordinate) const {
      return pUnderlyingHiRes->getIndexFromCoordinate(this->wrapCoordinate(coordinate - windowLowerCornerInclusive));
    }

    GridPtrType makeSubsampled(size_t ratio) const override {
      return std::make_shared<ResolutionMatchingGrid<T>>(pUnderlyingHiRes->makeSubsampled(ratio),
                                                         pUnderlyingLoRes->makeSubsampled(ratio));
    }

    GridPtrType makeSupersampled(size_t ratio) const override {
      return std::make_shared<ResolutionMatchingGrid<T>>(pUnderlyingHiRes->makeSupersampled(ratio),
                                                         pUnderlyingLoRes->makeSupersampled(ratio));
    }

    bool containsCellWithCoordinate(const Coordinate<int> &coord) const override {
      return this->pUnderlyingLoResInterpolated->containsCellWithCoordinate(coord);
    }

    bool containsCell(size_t i) const override {
      return this->pUnderlyingLoResInterpolated->containsCell(i);
    }

  };


  /*! \brief Virtual grid that points to a sub-section of another grid, rather than the whole grid.
    The virtual grid is interpreted as a subset of the underlying grid, at a position given by offset.
    Technically, there is no reason the virtual grid can't actually be larger than or
    lie partially outside the underlying grid, but
    if it is, flagged cells in the virtual grid won't propagate to the underlying grid unless they actually
    lie inside the underlying grid too.
  */
  template<typename T>
  class SectionOfGrid : public VirtualGrid<T> {
  private:
    Coordinate<int> cellOffset; //!< Co-ordinates of the lower left front corner of the sub-section (virtual grid) in the underlying grid.
    Coordinate<T> posOffset; //!< Position offset in co-moving co-ordinates of the lower front left corner of the sub-section.

  protected:
    using typename Grid<T>::GridPtrType;


  public:
    //! Converts an index of a cell in the virtual grid (subsection) into an index in the underlying grid.
    size_t mapIndexToUnderlying(size_t sec_id) const {
      auto coord = this->getCoordinateFromIndex(sec_id); // co-ordinate relative to virtual grid
      coord += cellOffset; // co-ordinate relative to underlying grid
      coord = this->wrapCoordinate(coord);
      if (!this->pUnderlying->containsCellWithCoordinate(coord))
        throw std::out_of_range("Out of range in SectionOfGrid::mapIndexToUnderlying");
      return this->pUnderlying->getIndexFromCoordinate(coord);
    }

    //! Converts an index of a cell in the *underlying* grid to a coordinate relative to *this* grid
    Coordinate<int> coordinateFromUnderlyingId(size_t underlying_id) const {
      Coordinate<int> coord = this->pUnderlying->getCoordinateFromIndex(
        underlying_id); // co-ordinate relative to underlying
      coord -= cellOffset; // co-ordinate relative to virtual grid
      coord = this->wrapCoordinate(coord);

      return coord;


    }

    /*! \brief Constructs a SectionOfGrid with the specified underlying grid, offset, and size

        \param pUnderlying - pointer to underlying grid.
        \param deltax - x component of a vector pointing from the underlying grid to the SectionOfGrid
        \param deltay - y component of a vector pointing from the underlying grid to the SectionOfGrid
        \param deltaz - z component of a vector pointing from the underlying grid to the SectionOfGrid

    */
    SectionOfGrid(GridPtrType pUnderlying, int deltax, int deltay, int deltaz, size_t size) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->periodicDomainSize,
                     size, // virtual grid may have a different size to the underlying grid
                     pUnderlying->cellSize,
                     pUnderlying->offsetLower.x +
                     deltax * pUnderlying->cellSize, // delta points from the underlying grid to the virtual grid
                     pUnderlying->offsetLower.y + deltay * pUnderlying->cellSize,
                     pUnderlying->offsetLower.z + deltaz * pUnderlying->cellSize,
                     pUnderlying->cellMassFrac, pUnderlying->cellSofteningScale),
      cellOffset(deltax, deltay, deltaz), // points from underlying to virtual grid
      posOffset(deltax * this->cellSize, deltay * this->cellSize, deltaz * this->cellSize) {

    }


    //! Returns true if the specified point lies in both grids
    bool containsPoint(const Coordinate<T> &coord) const override {
      return VirtualGrid<T>::containsPoint(coord) && this->pUnderlying->containsPoint(coord);
    }

    bool containsCellWithCoordinate(const Coordinate<int> &coord) const override {
      return VirtualGrid<T>::containsCellWithCoordinate(coord) &&
             this->pUnderlying->containsCellWithCoordinate(this->wrapCoordinate(coord + cellOffset));
    }

    bool containsCell(size_t i) const override {
      auto coord = this->getCoordinateFromIndex(i);
      return containsCellWithCoordinate(coord);
    }


    //! Outputs debug information
    void debugName(std::ostream &s) const override {
      s << "SectionOfGrid";
    }

    //! Stores the indices (relative to the virtual grid) of any flagged cells in the underlying grid that lie in both grids in targetArray
    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      for (size_t ptcl: underlyingArray) {
        auto coord = this->coordinateFromUnderlyingId(ptcl);
        if (this->containsCellWithCoordinate(coord))
          targetArray.push_back(this->getIndexFromCoordinate(coord));
      }
      tools::sortAndEraseDuplicate(targetArray);
    }

    //! Flags the specified cells (interpreted as indices in the virtual grid) in the underlying grid if they lie inside it.
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


  /*! \brief Virtual grid that has a lower resolution than the underlying grid.

  In contrast to the similar SuperSampleGrid, retrieved flags from the underlying grid are downscaled to the virtual grid, and cells flagged in the
  virtual grid are upscaled to the underlying grid.

  Field data is accessed by coarse graining the field on the underlying grid.
  */
  template<typename T>
  class SubSampleGrid : public VirtualGrid<T> {
  private:
    int factor; //!< Resolution is factor times lower than underlying grid
    int factor3; //!< Virtual grid contains factor3 times fewer points than the underlying grid

  protected:
    using typename Grid<T>::GridPtrType;

  public:
    /*! \brief Constructor for a SubSampleGrid

        \param pUnderlying - pointer to the underlying grid.
        \param factor - integer factor by which to reduce the resolution relative to the underlying grid.
    */
    SubSampleGrid(std::shared_ptr<Grid<T>> pUnderlying, int factor) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->periodicDomainSize, pUnderlying->size / factor,
                     pUnderlying->cellSize * factor, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                     pUnderlying->offsetLower.z, pUnderlying->cellMassFrac * pow(factor, 3.0),
                     pUnderlying->cellSofteningScale),
      factor(factor) {
      factor3 = factor * factor * factor;
    }

    //! Outputs debug information
    void debugName(std::ostream &s) const override {
      s << "SubSampleGrid";
    }

    //! Stores flags on the underlying grid in the target array, as indices defined on the virtual grid, by downscaling the underlying flag vector.
    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      std::vector<size_t> underlyingArray, downscaledArray;
      this->pUnderlying->getFlaggedCells(underlyingArray);
      Grid<T>::downscaleCellFlagVector(underlyingArray, downscaledArray, this->pUnderlying.get(), this);
      targetArray.insert(targetArray.end(), downscaledArray.begin(), downscaledArray.end());
    }

    //! Flags the specified cells in sourceArray (interpreted as indices on the virtual grid) by upscaling them to the underlying grid.
    void flagCells(const std::vector<size_t> &sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::upscaleCellFlagVector(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->flagCells(targetArray);
    }


    //! \brief Applies a specified operation to each cell in the underlying grid that corresponds to the specified cell in the virtual grid
    /*!
        \param id - cell in the virtual grid. Operation is applied to all cells in the underlying grid that are mapped to this by downscaling.
        \param callback - function to be applied to each of the sub-cells on the underlying grid.
    */
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

      int localFactor3 = 0;

      for (auto xi = coord0.x; xi < coord1.x; ++xi) {
        for (auto yi = coord0.y; yi < coord1.y; ++yi) {
          for (auto zi = coord0.z; zi < coord1.z; ++zi) {
            if(this->pUnderlying->containsCellWithCoordinate(Coordinate<int>(xi, yi, zi))) {
              callback(this->pUnderlying->getIndexFromCoordinateNoWrap(xi, yi, zi));
              localFactor3++;
            }
          }

        }
      }

      // At least one cell must be included, otherwise we are out of range:
      assert(localFactor3!=0);
      return localFactor3;
    }

    GridPtrType makeSupersampled(size_t ratio) const override {
      // Special case: supersampling a subsampled grid can needlessly destroy accuracy, so avoid it!
      if (ratio > factor) {
        // relative to the underlying grid we will be supersampled
        size_t ratioToUnderlying = tools::getRatioAndAssertInteger(static_cast<float>(ratio),
                                                                   static_cast<float>(factor));
        return this->pUnderlying->makeSupersampled(ratioToUnderlying);
      } else if (ratio == factor) {
        // gives back the underlying grid
        return std::const_pointer_cast<Grid<T>>(this->pUnderlying);
      } else {
        // relative to the underlying grid we will still be subsampled (just not as much, presumably)
        size_t ratioToUnderlying = tools::getRatioAndAssertInteger<float>(static_cast<float>(factor),
                                                                          static_cast<float>(ratio));
        return this->pUnderlying->makeSubsampled(ratioToUnderlying);
      }
    }

    GridPtrType makeSubsampled(size_t ratio) const override {
      return this->pUnderlying->makeSubsampled(ratio * factor);
    }

    virtual bool containsCellWithCoordinate(const Coordinate<int> &coord) const override {
      Coordinate<int> underlyingLowLeft = coord*factor;

      // Consider the cell to be included if at least one of the higher resolution cells is available:
      return this->pUnderlying->containsCellWithCoordinate(underlyingLowLeft) ||
             this->pUnderlying->containsCellWithCoordinate(underlyingLowLeft + (factor - 1));
    }

    bool containsCell(size_t i) const override {
      auto coord = this->getCoordinateFromIndex(i);
      return this->containsCellWithCoordinate(coord);
    }


  };


  /*! \brief Virtual grid which has a different mass fraction in each cell to the underlying grid.
    This is particularly useful for splitting the matter content into dark matter an baryons, as the
    cell mass is given by the total matter density parameter, and then two mass scaled grids are defined
    which describe the dark matter and baryons respectively, with the appropriate fraction of the total
    matter density.
  */
  template<typename T>
  class MassScaledGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::GridPtrType;


  public:
    //! \brief Constructor - specifies the underlying grid and the fraction to scale the mass by
    /*!
        \param pUnderlying - pointer to the underlying grid
        \param massScale - fraction to multiply the mass by
    */
    MassScaledGrid(GridPtrType pUnderlying, T massScale) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->periodicDomainSize, pUnderlying->size,
                     pUnderlying->cellSize, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                     pUnderlying->offsetLower.z, massScale * pUnderlying->cellMassFrac,
                     pUnderlying->cellSofteningScale) {}

    //! Outputs debug information
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
    const Coordinate<int> offset; //!< Vector; points from the centre of this virtual grid to the centre of the underlying grid.
  public:
    //! \brief Constructor - specified the underlying grid and the center point in co-moving co-ordinates corresponding to the centre
    /*!
        \param pUnderlying - pointer to the underlying grid
        \param center - point in co-moving co-ordinates corresponding to the centre
    */
    CenteredGrid(GridPtrType pUnderlying, Coordinate<T> center) :
      VirtualGrid<T>(pUnderlying,
                     pUnderlying->periodicDomainSize, pUnderlying->size,
                     pUnderlying->cellSize,
                     pUnderlying->offsetLower.x,
                     pUnderlying->offsetLower.y,
                     pUnderlying->offsetLower.z,
                     pUnderlying->cellMassFrac,
                     pUnderlying->cellSofteningScale),
      offset((Coordinate<T>(0.5 * this->pUnderlying->thisGridSize) - center) / this->pUnderlying->cellSize) {}

    //! Outputs debug information
    void debugName(std::ostream &s) const override {
      s << "CenteredGrid";
    }

    //! Returns the offset vector in co-moving co-ordinates.
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
    /*!
    Ie, points from the centre of the underlying grid to the centre of this grid, in
    integer units
    */
    Coordinate<int> getInverseOffset() const {
      return Coordinate<int>(-this->offset.x, -this->offset.y, -this->offset.z);
    }

    //! \brief Inverse offset, in Mpc/h
    Coordinate<T> getInversePointOffset() const {
      return this->getPointOffsetfromCoordinateOffset(this->getInverseOffset());
    }

  protected:

    //! \brief Returns index in underlying grid from co-ordinate relative to centred grid (unsigned int version)
    size_t getIndexFromCoordinateNoWrap(size_t x, size_t y, size_t z) const override {
      return this->pUnderlying->getIndexFromIndexAndStep(
        this->pUnderlying->getIndexFromCoordinateNoWrap(x, y, z), this->getInverseOffset());
    }

    //! \brief Returns index in underlying grid from co-ordinate relative to centred grid (int version)
    size_t getIndexFromCoordinateNoWrap(int x, int y, int z) const override {
      return this->pUnderlying->getIndexFromIndexAndStep(
        this->pUnderlying->getIndexFromCoordinateNoWrap(x, y, z), this->getInverseOffset());
    }

    //! \brief Returns index in underlying grid from co-ordinate relative to centred grid (coordinate version)
    size_t getIndexFromCoordinateNoWrap(const Coordinate<int> &coordinate) const override {
      return this->pUnderlying->getIndexFromIndexAndStep(
        this->pUnderlying->getIndexFromCoordinateNoWrap(coordinate), this->getInverseOffset());
    }

    //! \brief Returns index in underlying grid from co-ordinate relative to centred grid (without wrapping)
    size_t getIndexFromCoordinate(Coordinate<int> coord) const override {
      coord = this->pUnderlying->wrapCoordinate(coord);
      return this->getIndexFromCoordinateNoWrap(coord);
    }

    //! \brief Coordinate of a cell relative to the centred grid
    Coordinate<int> getCoordinateFromIndex(size_t id) const override {
      return this->pUnderlying->wrapCoordinate(this->pUnderlying->getCoordinateFromIndex(id) + this->offset);
    }

    //! \brief converts coordinates with respect to the underlying grid to co-ordinates with respect to the centred grid.
    Coordinate<T> getCentroidFromCoordinate(const Coordinate<int> &coord) const override {
      return this->pUnderlying->getCentroidFromCoordinate(
        this->pUnderlying->wrapCoordinate(this->offset + coord));
    }

    //! \brief Co-ordinate of centre of cell, with respect to centred grid
    Coordinate<T> getCentroidFromIndex(size_t id) const override {
      return this->wrapPoint(this->pUnderlying->getCentroidFromIndex(id) + this->getPointOffset());
    }

    //! \brief Returns index in underlying grid from co-ordinate relative to centred grid
    size_t getIndexFromPoint(Coordinate<T> point) const override {
      return this->pUnderlying->getIndexFromPoint(point + this->getInversePointOffset());
    }

    //! \brief Returns co-ordinate in underlying grid from co-ordinate relative to centred grid
    Coordinate<int> getCoordinateFromPoint(Coordinate<T> point) const override {
      return this->pUnderlying->getCoordinateFromPoint(point);
    }

    //! \brief Returns the flagged cells from the underlying grid.
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
    //! \brief Constructor - specifies underlying grid and the offset to shift by
    /*!
        \param pUnderlying - pointer to underlying grid
        \param dx - x component of offset
        \param dy - y component of offset
        \param dz - z component of offset
    */
    OffsetGrid(GridPtrType pUnderlying, T dx, T dy, T dz) :
      VirtualGrid<T>(pUnderlying, pUnderlying->periodicDomainSize, pUnderlying->size,
                     pUnderlying->cellSize,
                     pUnderlying->wrapIndividualCoordinate(pUnderlying->offsetLower.x + dx),
                     pUnderlying->wrapIndividualCoordinate(pUnderlying->offsetLower.y + dy),
                     pUnderlying->wrapIndividualCoordinate(pUnderlying->offsetLower.z + dz),
                     pUnderlying->cellMassFrac,
                     pUnderlying->cellSofteningScale) {}


    //! Outputs debug information
    void debugName(std::ostream &s) const override {
      s << "OffsetGrid";
    }
  };

  //! Wraps any VirtualGrid but allows it to have its own cell flags (independent of the original underlying grid)
  //! At construction, the new grid inherits the flags of the underlying grid, but these can then be manipualated
  //! fully independently.
  template<typename T>
  class IndependentFlaggingGrid : public VirtualGrid<T> {
  protected:
    using typename Grid<T>::GridPtrType;

  public:
    IndependentFlaggingGrid(const GridPtrType &pUnderlying) : VirtualGrid<T>(pUnderlying) {
      std::vector<size_t> ar;
      pUnderlying->getFlaggedCells(ar);
      this->flagCells(ar);
    }

    void flagCells(const std::vector<size_t> &sourceArray) override {
      Grid<T>::flagCells(sourceArray);
    }

    void unflagAllCells() override {
      Grid<T>::unflagAllCells();
    }

    size_t numFlaggedCells() const override {
      return Grid<T>::numFlaggedCells();
    }

    void getFlaggedCells(std::vector<size_t> &targetArray) const override {
      Grid<T>::getFlaggedCells(targetArray);
    }

    void debugInfo(std::ostream &s) const override {
      s << "[";
      this->pUnderlying->debugInfo(s);
      s << "]";
    }

    virtual GridPtrType withIndependentFlags() const override {
      return const_cast<IndependentFlaggingGrid<T> *>(this)->shared_from_this();
    }

    virtual GridPtrType withCoupledFlags() const override {
      return this->pUnderlying->withCoupledFlags();
    }

    virtual GridPtrType makeSupersampled(size_t ratio) const override {
      return std::make_shared<IndependentFlaggingGrid<T>>(this->pUnderlying->makeSupersampled(ratio));
    }

    virtual GridPtrType makeSubsampled(size_t ratio) const override {
      return std::make_shared<IndependentFlaggingGrid<T>>(this->pUnderlying->makeSubsampled(ratio));
    }
  };


}


#endif
