#ifndef IC_FIELD_HPP
#define IC_FIELD_HPP

#include <memory>
#include <vector>
#include <cassert>
#include <src/simulation/filters/filter.hpp>
#include <src/tools/numerics/fourier.hpp>
#include "src/io/numpy.hpp"
/*!
    \namespace fields
    \brief Define random fields on multiple grid levels
 */


#include "src/simulation/grid/grid.hpp"


// implementation in fourier.hpp
namespace tools {
  namespace numerics {
    namespace fourier {
      template<typename T>
      class FieldFourierManager;
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

  protected:
    const TPtrGrid pGrid;
    std::shared_ptr<FourierManager> fourierManager;
    TData data;
    bool fourier;

  public:
    //! Construct a field on the specified grid by moving the given field
    Field(Field<DataType, CoordinateType> &&move) : pGrid(move.pGrid), data(std::move(move.data)),
                                                    fourier(move.fourier) {
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size() == fourierManager->getRequiredDataSize());
    }

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

    TGrid &getGrid() const {
      return const_cast<TGrid &>(*pGrid);
    }

    TData &getDataVector() {
      return data;
    }

    const TData &getDataVector() const {
      return data;
    }

    operator std::vector<DataType> &() {
      return getDataVector();
    }

    operator const std::vector<DataType> &() const {
      return getDataVector();
    }


    DataType evaluateNearest(const Coordinate<CoordinateType> &location) const {
      auto offsetLower = pGrid->offsetLower;
      int x_p_0, y_p_0, z_p_0;


      // grid coordinates of parent cell that we're in
      x_p_0 = (int) floor(((location.x - offsetLower.x) / pGrid->cellSize));
      y_p_0 = (int) floor(((location.y - offsetLower.y) / pGrid->cellSize));
      z_p_0 = (int) floor(((location.z - offsetLower.z) / pGrid->cellSize));

      return (*this)[pGrid->getIndexFromCoordinate(Coordinate<int>(x_p_0, y_p_0, z_p_0))];
    }


    DataType evaluateInterpolated(Coordinate<CoordinateType> location) const {
      auto offsetLower = pGrid->offsetLower;
      int x_p_0, y_p_0, z_p_0, x_p_1, y_p_1, z_p_1;

      bool allowWrap = pGrid->coversFullSimulation();

      location -= offsetLower;
      location = pGrid->wrapPoint(location);

      // grid coordinates of parent cell starting to bottom-left
      // of our current point
      std::tie(x_p_0, y_p_0, z_p_0) = floor(location / pGrid->cellSize - 0.5);

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
    }

    const DataType &operator[](size_t i) const {
      return data[i];
    }

    DataType &operator[](size_t i) {
      return data[i];
    }

    /*
     * FOURIER TRANSFORM SUPPORT
     */

    template<typename T>
    void matchFourier(const T &other) {
      if (other.isFourier())
        toFourier();
      else
        toReal();
    }

    bool isFourier() const {
      return fourier;
    }

    void setFourier(bool fourier) {
      // Set a flag to indicate whether this field is in Fourier space or not - without actually applying any transform
      this->fourier = fourier;
    }

    void toFourier() {
      if (fourier) return;
      fourierManager->performTransform();
      assert(fourier);
    }

    void toReal() {
      if (!fourier) return;
      fourierManager->performTransform();
      assert(!fourier);
    }

    ComplexType getFourierCoefficient(int kx, int ky, int kz) const {
      return fourierManager->getFourierCoefficient(kx, ky, kz);
    }

    template<typename... Args>
    auto generateNewFourierFields(Args &&... args) {
      return fourierManager->generateNewFourierFields(args...);
    }

    //! Iterate (potentially in parallel) over each Fourier cell.
    /*!
     * The passed function takes arguments (value, kx, ky, kz) where value is the Fourier coeff value
     * at k-mode kx, ky, kz.
     * If the function returns a value, the Fourier coeff in that cell is updated accordingly.
     */
    template<typename... Args>
    void forEachFourierCell(Args &&... args) {
      fourierManager->forEachFourierCell(args...);
    }

    template<typename... Args>
    void forEachFourierCell(Args &&... args) const {
      fourierManager->forEachFourierCell(args...);
    }

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

    void applyFilter(const filters::Filter<CoordinateType> &filter) {
      forEachFourierCell([&filter](ComplexType current_value, CoordinateType kx, CoordinateType ky, CoordinateType kz) {
        CoordinateType k = sqrt(double(kx * kx + ky * ky + kz * kz));
        return current_value * filter(k);
      });
    }


    void addFieldFromDifferentGrid(const Field<DataType, CoordinateType> &source) {
      assert(!source.isFourier());
      toReal();
      TPtrGrid pSourceProxyGrid = source.getGrid().makeProxyGridToMatch(getGrid());

      auto evaluator = makeEvaluator(source, *pSourceProxyGrid);
      evaluator->addTo(*this);

    }


#ifdef FILTER_ON_COARSE_GRID
    // Version for compatibility with pre-Dec 2016 output
    // Applies filter BEFORE interpolating onto the fine grid, which results in more pixel window function
    // artefacts but (presumably?) fewer artefacts from the grid-level window function

    void addFieldFromDifferentGridWithFilter(Field<DataType, CoordinateType> & source,
                                             const filters::Filter<CoordinateType> & filter) {


      Field<DataType, CoordinateType> temporaryField(source);
      temporaryField.applyFilter(filter);

      auto & temporaryFieldData=temporaryField.data;

      size_t size3 = getGrid().size3;

      temporaryField.toReal();
      this->toReal();

#pragma omp parallel for schedule(static)
      for (size_t ind_l = 0; ind_l < size3; ind_l++) {
        data[ind_l] += temporaryField.evaluateInterpolated(pGrid->getPointFromId(ind_l));
      }

    }

#else

    void addFieldFromDifferentGridWithFilter(Field<DataType, CoordinateType> &source,
                                             const filters::Filter<CoordinateType> &filter) {

      source.toReal();
      Field<DataType, CoordinateType> temporaryField(getGrid(), false);
      auto &temporaryFieldData = temporaryField.getDataVector();

      assert(data.size() == temporaryFieldData.size());

      size_t size = pGrid->size3;

#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < size; ++i) {
        temporaryFieldData[i] = source.evaluateInterpolated(pGrid->getCentroidFromIndex(i));
      }

      temporaryField.applyFilter(filter);

      this->matchFourier(temporaryField); // expect that the temporary field is now stored in Fourier space

#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < temporaryFieldData.size(); ++i) {
        data[i] += temporaryFieldData[i];
      }
    }

#endif

    void addFieldFromDifferentGrid(Field<DataType, CoordinateType> &source) {
      source.toReal();
      addFieldFromDifferentGrid(const_cast<const Field<DataType, CoordinateType> &>(source));
    }

    void dumpGridData(std::string filename) const {
      int n = static_cast<int>(getGrid().size);
      const int dim[3] = {n, n, n};
      io::numpy::SaveArrayAsNumpy(filename, false, 3, dim, data.data());
    }

    void loadGridData(std::string filename) {
      int n = static_cast<int>(getGrid().size);
      int n0, n1, n2;
      io::numpy::LoadArrayFromNumpy(filename, n0, n1, n2, data );
      if(n0!=n || n1!=n || n2!=n) {
        throw std::runtime_error("Incorrect size for imported numpy array");
      }
      assert(data.size()==getGrid().size3);
      data.resize(fourierManager->getRequiredDataSize());
    }


  };

  //! Mostly used for debugging. Convert a field from holding one data type to another, e.g. complex to real.
  template<typename TargetDataType, typename SourceDataType, typename CoordinateType>
  std::shared_ptr<Field<TargetDataType, CoordinateType>>
  convertField(Field<SourceDataType, CoordinateType> &field) {
    assert(!field.isFourier());
    auto newField = make_shared<Field<TargetDataType, CoordinateType>>(field.getGrid(), field.isFourier());
    auto &newData = newField->getDataVector();
    auto &originalData = field.getDataVector();
    size_t size = field.getGrid().size3;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < size; ++i) {
      newData[i] = tools::datatypes::real_part_if_complex(originalData[i]);
    }

    return newField;
  };

  template<typename TargetDataType, typename CoordinateType>
  std::shared_ptr<Field<TargetDataType, CoordinateType>>
  convertField(Field<TargetDataType, CoordinateType> &field) {
    return field.shared_from_this();
  };


}


#endif
