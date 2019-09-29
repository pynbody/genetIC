#ifndef IC_FIELD_HPP
#define IC_FIELD_HPP

#include <memory>
#include <vector>
#include <cassert>
#include <src/simulation/filters/filter.hpp>
#include <src/tools/numerics/fourier.hpp>
#include "src/io/numpy.hpp"
#include "src/simulation/grid/grid.hpp"
#include "src/tools/numerics/tricubic.hpp"


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
      template<typename T>
      class FieldFourierManager;
      // implementation in fourier.hpp
      
    }
  }
}


namespace fields {

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

  //! Class to manage and evaluate a field defined on a single grid.
  template<typename DataType, typename CoordinateType=tools::datatypes::strip_complex<DataType>>
  class Field : public std::enable_shared_from_this<Field<DataType, CoordinateType>> {
  public:
    using TGrid = const grids::Grid<CoordinateType>;
    using TPtrGrid = std::shared_ptr<TGrid>;
    using TData = std::vector<DataType>;
    using value_type = DataType;
    using ComplexType = tools::datatypes::ensure_complex<DataType>;

    using FourierManager = tools::numerics::fourier::FieldFourierManager<DataType>;
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
    Field(Field<DataType, CoordinateType> &&move) : pGrid(move.pGrid), data(std::move(move.data)),
                                                    fourier(move.fourier) {
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());
    }

    //! Copy constructor
    Field(const Field<DataType, CoordinateType> &copy)
      : std::enable_shared_from_this<Field<DataType, CoordinateType>>(),
        pGrid(copy.pGrid), data(copy.data),
        fourier(copy.fourier) {
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());
    }

    //! Construct a field on the specified grid by moving the given data
    Field(TGrid &grid, TData &&dataVector, bool fourier = true) : pGrid(grid.shared_from_this()),
                                                                  data(std::move(dataVector)), fourier(fourier) {

      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());

    }

    //! Construct a field on the specified grid by copying the given data
    Field(TGrid &grid, const TData &dataVector, bool fourier = true) : pGrid(grid.shared_from_this()),
                                                                       data(dataVector), fourier(fourier) {

      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());

    }

    //! Construct a zero-filled field on the specified grid
    Field(TGrid &grid, bool fourier = true) : pGrid(grid.shared_from_this()),
                                              fourierManager(std::make_shared<FourierManager>(*this)),
                                              data(fourierManager->getRequiredDataSize(), 0),
                                              fourier(fourier) {
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

    //! Add the specified value to the field and the given location, using conjugate deinterpolation
    /*! For an explanation of what is meant by 'conjugate deinterpolation' see the
     * documentation for numerics::LocalUnitTricubicApproximation::getTransposeElementsForPosition
     */
    void deInterpolate(Coordinate<CoordinateType> location, DataType value) {
      int x_p_0, y_p_0, z_p_0;
      CoordinateType dx, dy, dz;

      location -= pGrid->offsetLower;
      location = pGrid->wrapPoint(location);

#ifdef CUBIC_INTERPOLATION
      // grid coordinates of parent cell whose *centroid* (not corner) is to the bottom-left of our current point
      std::tie(x_p_0, y_p_0, z_p_0) = floor(location / pGrid->cellSize - 0.5);

      DataType valsForInterpolation[4][4][4];
      std::tie(dx,dy,dz) = (location / pGrid->cellSize - 0.5);
      dx-=x_p_0;
      dy-=y_p_0;
      dz-=z_p_0;

      // Caching would hugely speed this up, since we anticipate repeated calls with the same dx,dy,dz (to numerical accuracy)
    numerics::LocalUnitTricubicApproximation<DataType>::getTransposeElementsForPosition(dx,dy,dz,valsForInterpolation);

      for(int i=-1; i<3; ++i) {
        for(int j=-1; j<3; ++j) {
          for(int k=-1; k<3; ++k) {
            (*this)[pGrid->getIndexFromCoordinate({x_p_0+i, y_p_0+j, z_p_0+k})] += value * valsForInterpolation[i+1][j+1][k+1];
          }
        }
      }
#else 
    // TODO: The below is not conjugate deinterpolation for the linear interpolation
    // implemented when CUBIC_INTERPOLATION is off. However, switching off 
    // CUBIC_INTERPOLATION is not recommended, so fixing this inconsistency is low
    // priority. In reality, the deinterpolation here is conjugate to zero-order
    // (nearest neighbour) interpolation.
    std::tie(x_p_0, y_p_0, z_p_0) = floor(location / pGrid->cellSize);
    (*this)[pGrid->getIndexFromCoordinate({x_p_0, y_p_0, z_p_0})] += value;    
#endif
    }

    //! Evaluates the field at the specified co-ordinate using interpolation.
    DataType evaluateInterpolated(Coordinate<CoordinateType> location) const {
      int x_p_0, y_p_0, z_p_0;

      location -= pGrid->offsetLower;
      location = pGrid->wrapPoint(location);

      // grid coordinates of parent cell whose *centroid* (not corner) is to the bottom-left of our current point
      std::tie(x_p_0, y_p_0, z_p_0) = floor(location / pGrid->cellSize - 0.5);

#ifdef CUBIC_INTERPOLATION
      // The following cubic interpolation constructor is potentially expensive. If it turns out to be a major
      // part of the overall runtime, a caching scheme could be implemented to reduce the overall cost of
      // interpolation, especially when interpolating an entire grid.
      DataType valsForInterpolation[4][4][4];
      for(int i=-1; i<3; ++i) {
        for(int j=-1; j<3; ++j) {
          for(int k=-1; k<3; ++k) {
            valsForInterpolation[i+1][j+1][k+1] = (*this)[pGrid->getIndexFromCoordinate({x_p_0+i, y_p_0+j, z_p_0+k})];
          }
        }
      }
      numerics::LocalUnitTricubicApproximation<DataType> interpolator(valsForInterpolation);
      CoordinateType dx,dy,dz;

      // Work out the fractional displacement of our target point between the centroid of the cell identified above
      // and the next one along. dx, dy, dz will be between zero and one unless something goes badly wrong!
      std::tie(dx,dy,dz) = (location / pGrid->cellSize - 0.5);
      dx-=x_p_0;
      dy-=y_p_0;
      dz-=z_p_0;

      return interpolator(dx,dy,dz);
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

#pragma omp parallel for schedule(static) default(none)  shared(data)
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
