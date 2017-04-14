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

  /** Class to manage and evaluate a field defined on a single grid.  */
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
    Field(Field<DataType, CoordinateType> &&move) : pGrid(move.pGrid), data(std::move(move.data)),
                                                    fourier(move.fourier) {
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size()==fourierManager->getRequiredDataSize());

    }

    Field(TGrid &grid, TData &&dataVector, bool fourier = true) : pGrid(grid.shared_from_this()),
                                                                  data(std::move(dataVector)), fourier(fourier) {
      /*
       * Construct a field on the specified grid by moving the given data
       */
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size()==fourierManager->getRequiredDataSize());

    }

    Field(TGrid &grid, const TData &dataVector, bool fourier = true) : pGrid(grid.shared_from_this()),
                                                                       data(dataVector), fourier(fourier) {
      /*
       * Construct a field on the specified grid by copying the given data
       */
      fourierManager = std::make_shared<FourierManager>(*this);
      assert(data.size()==fourierManager->getRequiredDataSize());

    }

    Field(TGrid &grid, bool fourier = true) : pGrid(grid.shared_from_this()),
                                              fourierManager(std::make_shared<FourierManager>(*this)),
                                              data(fourierManager->getRequiredDataSize(), 0),
                                              fourier(fourier) {
      /*
       * Construct a zero-filled field on the specified grid
       */
    }

  public:

    TGrid &getGrid() const {
      return const_cast<TGrid &>(*pGrid);
    }

    /*
     * ACCESS TO UNDERLYING DATA VECTOR
     */

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

    /*
     * EVALUATION OPERATIONS
     */

    DataType evaluateNearest(const Coordinate<CoordinateType> &location) const {
      auto offsetLower = pGrid->offsetLower;
      int x_p_0, y_p_0, z_p_0;


      // grid coordinates of parent cell that we're in
      x_p_0 = (int) floor(((location.x - offsetLower.x) / pGrid->dx));
      y_p_0 = (int) floor(((location.y - offsetLower.y) / pGrid->dx));
      z_p_0 = (int) floor(((location.z - offsetLower.z) / pGrid->dx));

      return (*this)[pGrid->getCellIndex(Coordinate<int>(x_p_0, y_p_0, z_p_0))];
    }


    DataType evaluateInterpolated(const Coordinate<CoordinateType> &location) const {
      auto offsetLower = pGrid->offsetLower;
      int x_p_0, y_p_0, z_p_0, x_p_1, y_p_1, z_p_1;


      // grid coordinates of parent cell starting to bottom-left
      // of our current point
      x_p_0 = (int) floor(((location.x - offsetLower.x) / pGrid->dx - 0.5));
      y_p_0 = (int) floor(((location.y - offsetLower.y) / pGrid->dx - 0.5));
      z_p_0 = (int) floor(((location.z - offsetLower.z) / pGrid->dx - 0.5));

      // grid coordinates of top-right
      x_p_1 = x_p_0 + 1;
      y_p_1 = y_p_0 + 1;
      z_p_1 = z_p_0 + 1;

      // weights, which are the distance to the centre point of the
      // upper-right cell, in grid units (-> maximum 1)

      CoordinateType xw0, yw0, zw0, xw1, yw1, zw1;
      xw0 = ((CoordinateType) x_p_1 + 0.5) - ((location.x - offsetLower.x) / pGrid->dx);
      yw0 = ((CoordinateType) y_p_1 + 0.5) - ((location.y - offsetLower.y) / pGrid->dx);
      zw0 = ((CoordinateType) z_p_1 + 0.5) - ((location.z - offsetLower.z) / pGrid->dx);

      xw1 = 1. - xw0;
      yw1 = 1. - yw0;
      zw1 = 1. - zw0;

      assert(xw0 <= 1.0 && xw0 >= 0.0);

      // allow things on the boundary to 'saturate' value, but beyond boundary
      // is not acceptable
      //
      // TODO - in some circumstances we may wish to replace this with wrapping
      // but not all circumstances!
      int size_i = static_cast<int>(pGrid->size);
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

      // return (*this)[pGrid->getCellIndexNoWrap(x_p_0, y_p_0, z_p_0)];


      return xw0 * yw0 * zw1 * (*this)[pGrid->getCellIndexNoWrap(x_p_0, y_p_0, z_p_1)] +
             xw1 * yw0 * zw1 * (*this)[pGrid->getCellIndexNoWrap(x_p_1, y_p_0, z_p_1)] +
             xw0 * yw1 * zw1 * (*this)[pGrid->getCellIndexNoWrap(x_p_0, y_p_1, z_p_1)] +
             xw1 * yw1 * zw1 * (*this)[pGrid->getCellIndexNoWrap(x_p_1, y_p_1, z_p_1)] +
             xw0 * yw0 * zw0 * (*this)[pGrid->getCellIndexNoWrap(x_p_0, y_p_0, z_p_0)] +
             xw1 * yw0 * zw0 * (*this)[pGrid->getCellIndexNoWrap(x_p_1, y_p_0, z_p_0)] +
             xw0 * yw1 * zw0 * (*this)[pGrid->getCellIndexNoWrap(x_p_0, y_p_1, z_p_0)] +
             xw1 * yw1 * zw0 * (*this)[pGrid->getCellIndexNoWrap(x_p_1, y_p_1, z_p_0)];
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

    void setFourierCoefficient(int kx, int ky, int kz, const ComplexType & value) {
      fourierManager->setFourierCoefficient(kx, ky, kz, value);
    }

    void applyFilter(const filters::Filter<CoordinateType> &filter) {
      toFourier();
      size_t size3 = pGrid->size3;
      int size = static_cast<int>(pGrid->size);


      CoordinateType kMin = getGrid().getFourierKmin();
      for(int kx=0; kx<size/2+1; kx++) {
        int ky_lower = -size/2+1;
        if(kx==0) ky_lower = 0;
        for(int ky=ky_lower; ky<size/2+1; ky++) {
          int kz_lower = -size/2+1;
          if(kx==0 && ky==0)
            kz_lower = 0;

          for(int kz=kz_lower; kz<size/2+1; kz++) {
            auto coeff = getFourierCoefficient(kx, ky, kz);
            CoordinateType k = kMin*sqrt(double(kx*kx+ky*ky+kz*kz));
            coeff*=filter(k);
            setFourierCoefficient(kx,ky,kz, coeff);
          }
        }
      }
    }


    void addFieldFromDifferentGrid(const Field<DataType, CoordinateType> &source) {
      assert(!source.isFourier());
      toReal();
      TPtrGrid pSourceProxyGrid = source.getGrid().makeProxyGridToMatch(getGrid());

      size_t size3 = getGrid().size3;

#pragma omp parallel for schedule(static)
      for (size_t ind_l = 0; ind_l < size3; ind_l++) {
        if (pSourceProxyGrid->containsCell(ind_l))
          data[ind_l] += pSourceProxyGrid->getFieldAt(ind_l, source);
      }


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
        data[ind_l] += temporaryField.evaluateInterpolated(pGrid->getCellCentroid(ind_l));
      }

    }

#else

    void addFieldFromDifferentGridWithFilter(Field<DataType, CoordinateType> &source,
                                             const filters::Filter<CoordinateType> &filter) {

      source.toReal();
      Field<DataType, CoordinateType> temporaryField(getGrid(), false);
      auto &temporaryFieldData = temporaryField.getDataVector();

      assert(data.size() == temporaryFieldData.size());

      size_t size = source.getGrid().size3;

#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < size; ++i) {
        temporaryFieldData[i] = source.evaluateInterpolated(pGrid->getCellCentroid(i));
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


  };


  template<typename CoordinateType>
  std::shared_ptr<Field<CoordinateType, CoordinateType>> getRealPart(Field<CoordinateType, CoordinateType> &field) {
    return field.shared_from_this();
  };


  template<typename CoordinateType>
  std::shared_ptr<Field<CoordinateType, CoordinateType>>
  getRealPart(Field<std::complex<CoordinateType>, CoordinateType> &field) {
    assert(!field.isFourier());
    auto realPart = make_shared<Field<CoordinateType, CoordinateType>>(field.getGrid(), field.isFourier());
    auto &realPartData = realPart->getDataVector();
    auto &complexData = field.getDataVector();
    size_t size = field.getGrid().size3;

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < size; ++i) {
      realPartData[i] = complexData[i].real();
    }

    return realPart;


  };


}


#endif //IC_FIELD_HPP
