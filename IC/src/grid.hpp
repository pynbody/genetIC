#ifndef __GRID_HPP
#define __GRID_HPP

#include <cassert>
#include <set>
#include <type_traits>
#include <memory>
#include <vector>
#include <complex>
#include "fft.hpp"
#include "filter.hpp"
#include "coordinate.hpp"

using namespace std;

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

  using TField = std::vector<std::complex<T>>;
  using TRealField = std::vector<T>;
  using PtrTField = std::shared_ptr<std::vector<std::complex<T>>>;
  using GridPtrType = std::shared_ptr<Grid<T>>;
  using ConstGridPtrType = std::shared_ptr<const Grid<T>>;

private:
  std::shared_ptr<std::vector<std::complex<T>>> pField;
  bool fieldFourier; //< is the field in k-space (true) or x-space (false)?

  // The grid offsets after Zeldovich approximation is applied
  // (nullptr before that):
  std::shared_ptr<TRealField> pOff_x;
  std::shared_ptr<TRealField> pOff_y;
  std::shared_ptr<TRealField> pOff_z;

  T hFactor, cellMass;
  T kMin, kMinSquared;

  std::vector<size_t> particleArray; // just a list of particles on this grid for one purpose or another
  std::vector<T *> particleProperties; // a list of particle properties

protected:

  static void upscaleParticleList(const std::vector<size_t> sourceArray,
                                  std::vector<size_t> &targetArray,
                                  const Grid<T> *source,
                                  const Grid<T> *target) {
    int x0, y0, z0, x1, y1, z1, x, y, z;

    assert(target->size >= source->size);
    assert((target->size) % (source->size) == 0);
    int factor = int(target->size / source->size);
    targetArray.clear();

    for (auto id: sourceArray) {
      auto coord = source->getCellCoordinate(id);
      iterateOverCube<int>(
        coord * factor, coord * factor + factor,
        [&targetArray, &target](const Coordinate<int> & subCoord) {
          targetArray.push_back(target->getCellIndexNoWrap(subCoord));
        }
      );
    }
  }

  static void downscaleParticleList(const std::vector<size_t> sourceArray,
                                    std::vector<size_t> &targetArray,
                                    const Grid<T> *source,
                                    const Grid<T> *target) {


    std::set<size_t> targetSet;
    int x, y, z;

    assert(source->size >= target->size);
    assert((source->size) % (target->size) == 0);
    size_t factor = source->size / target->size;

    for (auto id: sourceArray) {
      auto coord = source->getCellCoordinate(id);
      targetSet.insert(target->getCellIndexNoWrap(coord / factor));
    }
    targetArray.clear();
    targetArray.insert(targetArray.end(), targetSet.begin(), targetSet.end());
  }

  void setKmin() {
    kMin = 2. * M_PI / boxsize;
    kMinSquared = kMin * kMin;
  }

  static int getRatioAndAssertInteger(T p, T q) {
    const T tolerance = 1e-6;
    T ratio = p / q;
    int rounded_ratio = int(round(ratio));
    assert(abs(T(rounded_ratio) - ratio) < tolerance);
    return rounded_ratio;
  }

public:

  const T simsize, boxsize, dx;
  const Coordinate<T> offsetLower;
  const size_t size, size2, size3;


  Grid(T simsize, size_t n, T dx = 1.0, T x0 = 0.0, T y0 = 0.0, T z0 = 0.0, bool withField = true) :
    simsize(simsize), boxsize(dx * n),
    dx(dx), offsetLower(x0, y0, z0),
    size(n), size2(n * n), size3(n * n * n) {
    // cerr << "Grid ctor " << this <<  endl;
    if (withField) {
      pField = std::make_shared<std::vector<std::complex<T>>>(size3, 0);
      pField->shrink_to_fit();
    }
    setKmin();
  }


  Grid(size_t n) : simsize(0), boxsize(n),
                   dx(1.0), offsetLower(0, 0, 0),
                   size(n), size2(n * n), size3(n * n * n) {
    // cerr << "Grid ctor size-only" << endl;
    pField = nullptr;
    setKmin();
  }

  virtual ~Grid() {
    // cerr << "~Grid " << this << endl;
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


  GridPtrType makeProxyGridToMatch(const Grid<T> &target) const {
    GridPtrType proxy = std::const_pointer_cast<Grid<T>>(this->shared_from_this());
    if (target.dx > dx) {
      int ratio = getRatioAndAssertInteger(target.dx, dx);
      proxy = std::make_shared<SubSampleGrid<T>>(proxy, ratio);
    } else if (target.dx < dx) {
      int ratio = getRatioAndAssertInteger(dx, target.dx);
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

  template<typename TArray>
  void addFieldFromDifferentGrid(const Grid<T> &grid_src, const TArray &pField_x_src, TArray &pField_x_dest) {

    // TODO: make signature above const-correct

    GridPtrType pSourceProxyGrid = grid_src.makeProxyGridToMatch(*this);

    assert(pSourceProxyGrid->fieldIsSuitableSize(pField_x_src));
    assert(fieldIsSuitableSize(pField_x_dest));
    assert(pSourceProxyGrid->size3 == size3);


#pragma omp parallel for schedule(static)
    for (size_t ind_l = 0; ind_l < size3; ind_l++) {
      if (pSourceProxyGrid->containsCell(ind_l))
        pField_x_dest[ind_l] += pSourceProxyGrid->getFieldAt(ind_l, pField_x_src);

    }
  }

  void addFieldFromDifferentGrid(const Grid<T> &grid_src) {
    TField &pField_dest = getFieldReal();
    assert(!grid_src.isFieldFourier());
    const TField &pField_src = grid_src.getField();

    addFieldFromDifferentGrid(grid_src, pField_src, pField_dest);

    auto offset_fields_src = grid_src.getOffsetFields();
    auto offset_fields_dest = getOffsetFields();

    if (std::get<0>(offset_fields_src)->size() > 0 && std::get<0>(offset_fields_dest)->size() > 0) {
      addFieldFromDifferentGrid(grid_src, *std::get<0>(offset_fields_src), *std::get<0>(offset_fields_dest));
      addFieldFromDifferentGrid(grid_src, *std::get<1>(offset_fields_src), *std::get<1>(offset_fields_dest));
      addFieldFromDifferentGrid(grid_src, *std::get<2>(offset_fields_src), *std::get<2>(offset_fields_dest));
    }
  }

  void addFieldFromDifferentGrid(Grid<T> &grid_src) {
    grid_src.getFieldReal();
    addFieldFromDifferentGrid(const_cast<const Grid<T> &>(grid_src));
  }

  virtual void debugInfo(std::ostream &s) const {
    s << "Grid of side " << size << " address " << this << "; " << this->particleArray.size() << " cells marked";
  }


  virtual void gatherParticleList(std::vector<size_t> &targetArray) const {
    targetArray.insert(targetArray.end(), particleArray.begin(), particleArray.end());
  }

  virtual void distributeParticleList(const std::vector<size_t> &sourceArray) {
    particleArray.insert(particleArray.end(), sourceArray.begin(), sourceArray.end());
  }

  virtual void clearParticleList() {
    particleArray.clear();
  }

  virtual size_t estimateParticleListSize() {
    return particleArray.size();
  }


  virtual bool pointsToGrid(Grid<T> *pOther) {
    return this == pOther;
  }

  bool pointsToAnyGrid(std::vector<std::shared_ptr<Grid<T>>> grids) {
    for (auto g: grids) {
      if (pointsToGrid(g.get()))
        return true;
    }
    return false;
  }

  ///////////////////////////////////////
  //  Field manipulation routines
  ///////////////////////////////////////

  virtual TField &getFieldFourier() {
    assert(pField != nullptr);
    if (!fieldFourier) {
      fft(pField->data(), pField->data(), size, 1);
      fieldFourier = true;
    }
    return *pField;
  }

  virtual TField &getFieldReal() {
    assert(pField != nullptr);
    if (fieldFourier) {
      fft(pField->data(), pField->data(), size, -1);
      fieldFourier = false;
    }
    return *pField;
  }

  virtual const TField &getField() const {
    assert(pField != nullptr);
    return *pField;
  }

  TField &getField() {
    return const_cast<TField &>(const_cast<const Grid *>(this)->getField());
  }

  void applyFilter(const Filter<T> &filter, TField &fieldFourier) {
    assert(fieldFourier.size() == size3);

#pragma omp parallel for
    for (size_t i = 0; i < size3; ++i) {
      fieldFourier[i] *= filter(getFourierCellAbsK(i));
    }
  }


  virtual bool fieldIsSuitableSize(const TField &field) {
    return field.size() == size3;
  }

  virtual bool fieldIsSuitableSize(const TRealField &field) {
    return field.size() == size3;
  }

  virtual bool containsCell(const Coordinate<T> &coord) const {
    return coord.x >= 0 && coord.y >= 0 && coord.z >= 0 &&
           coord.x < size && coord.y < size && coord.z < size;
  }

  virtual bool containsCell(size_t i) const {
    return i<size3;
  }

  virtual complex<T> getFieldAt(size_t i, const TField &field) {
    return field[i];
  }

  virtual T getFieldAt(size_t i, const TRealField &field) {
    return field[i];
  }

  virtual complex<T> getFieldAt(size_t i) {
    return getFieldAt(i, getFieldReal());
  }


  auto getOffsetFields() {
    return std::make_tuple(pOff_x, pOff_y, pOff_z);
  }

  auto getOffsetFields() const {
    return std::make_tuple(pOff_x, pOff_y, pOff_z);
  }

  virtual bool isFieldFourier() const {
    return fieldFourier;
  }

  void getParticle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    getParticleNoWrap(id, x, y, z, vx, vy, vz, cellMassi, eps);
    simWrap(x, y, z);
  }

  virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    T x0, y0, z0;
    getParticleNoOffset(id, x, y, z, vx, vy, vz, cellMassi, eps);
    std::tie(x0, y0, z0) = getCellCentroid(id);
    x += x0;
    y += y0;
    z += z0;
  }

  virtual T getMass() const {
    return cellMass;
  }

  virtual T getEps() const {
    return dx * 0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
  }

  void getParticleNoOffset(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {

    x = (*pOff_x)[id];
    y = (*pOff_y)[id];
    z = (*pOff_z)[id];

    vx = (*pOff_x)[id] * hFactor;
    vy = (*pOff_y)[id] * hFactor;
    vz = (*pOff_z)[id] * hFactor;

    cellMassi = getMass();
    eps = getEps();
  }

  virtual void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {

    vx = getFieldInterpolated(x, y, z, *pOff_x);
    vy = getFieldInterpolated(x, y, z, *pOff_y);
    vz = getFieldInterpolated(x, y, z, *pOff_z);
    x += vx;
    y += vy;
    z += vz;

    simWrap(x, y, z);

    vx *= hFactor;
    vy *= hFactor;
    vz *= hFactor;

    cellMassi = getMass();
    eps = getEps();
  }


  template<typename TField>
  typename TField::value_type getFieldInterpolated(const Coordinate<T> &location, const TField &pField) const {
    return getFieldInterpolated(location.x, location.y, location.z, pField);
  }

  template<typename TField>
  typename TField::value_type getFieldInterpolated(const T &x, const T &y, const T &z, const TField &pField) const {
    int x_p_0, y_p_0, z_p_0, x_p_1, y_p_1, z_p_1;

    // grid coordinates of parent cell starting to bottom-left
    // of our current point
    x_p_0 = (int) floor(((x - offsetLower.x) / dx - 0.5));
    y_p_0 = (int) floor(((y - offsetLower.y) / dx - 0.5));
    z_p_0 = (int) floor(((z - offsetLower.z) / dx - 0.5));

    // grid coordinates of top-right
    x_p_1 = x_p_0 + 1;
    y_p_1 = y_p_0 + 1;
    z_p_1 = z_p_0 + 1;

    // weights, which are the distance to the centre point of the
    // upper-right cell, in grid units (-> maximum 1)

    T xw0, yw0, zw0, xw1, yw1, zw1;
    xw0 = ((T) x_p_1 + 0.5) - ((x - offsetLower.x) / dx);
    yw0 = ((T) y_p_1 + 0.5) - ((y - offsetLower.y) / dx);
    zw0 = ((T) z_p_1 + 0.5) - ((z - offsetLower.z) / dx);

    xw1 = 1. - xw0;
    yw1 = 1. - yw0;
    zw1 = 1. - zw0;

    assert(xw0 <= 1.0 && xw0 >= 0.0);

    // allow things on the boundary to 'saturate' value, but beyond boundary
    // is not acceptable
    //
    // TODO - in some circumstances we may wish to replace this with wrapping
    // but not all circumstances!
    int size_i = static_cast<int>(size);
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


    return xw0 * yw0 * zw1 * pField[getCellIndexNoWrap(x_p_0, y_p_0, z_p_1)] +
           xw1 * yw0 * zw1 * pField[getCellIndexNoWrap(x_p_1, y_p_0, z_p_1)] +
           xw0 * yw1 * zw1 * pField[getCellIndexNoWrap(x_p_0, y_p_1, z_p_1)] +
           xw1 * yw1 * zw1 * pField[getCellIndexNoWrap(x_p_1, y_p_1, z_p_1)] +
           xw0 * yw0 * zw0 * pField[getCellIndexNoWrap(x_p_0, y_p_0, z_p_0)] +
           xw1 * yw0 * zw0 * pField[getCellIndexNoWrap(x_p_1, y_p_0, z_p_0)] +
           xw0 * yw1 * zw0 * pField[getCellIndexNoWrap(x_p_0, y_p_1, z_p_0)] +
           xw1 * yw1 * zw0 * pField[getCellIndexNoWrap(x_p_1, y_p_1, z_p_0)];
  }


  virtual std::shared_ptr<Grid<T>> makeScaledMassVersion(T massRatio) {
    return std::make_shared<MassScaledGrid<T>>(this->shared_from_this(), massRatio);
  }


  virtual void zeldovich(T hfac, T particlecellMass) {

    hFactor = hfac;
    cellMass = particlecellMass;

    cout << "Applying Zeldovich approximation; grid cell size=" << dx << " Mpc/h...";
    cout.flush();

    // make three arrays for manipulating in fourier space
    auto psift1k = std::vector<std::complex<T>>(size3, 0);
    auto psift2k = std::vector<std::complex<T>>(size3, 0);
    auto psift3k = std::vector<std::complex<T>>(size3, 0);

    // get a reference to the density field in fourier space
    auto &pField_k = getFieldFourier();

    int iix, iiy, iiz;
    T kfft;
    size_t idx;

    T kw = 2. * M_PI / boxsize;

#pragma omp parallel for schedule(static) default(shared) private(iix, iiy, iiz, kfft, idx)
    for (size_t ix = 0; ix < size; ix++) {
      for (size_t iy = 0; iy < size; iy++) {
        for (size_t iz = 0; iz < size; iz++) {

          idx = (ix * size + iy) * size + iz;

          if (ix > size / 2) iix = static_cast<int>(ix) - size; else iix = ix;
          if (iy > size / 2) iiy = static_cast<int>(iy) - size; else iiy = iy;
          if (iz > size / 2) iiz = static_cast<int>(iz) - size; else iiz = iz;

          kfft = (T) (iix * iix + iiy * iiy + iiz * iiz);

          psift1k[idx].real(-pField_k[idx].imag() / (T) (kfft) * iix / kw);
          psift1k[idx].imag(pField_k[idx].real() / (T) (kfft) * iix / kw);
          psift2k[idx].real(-pField_k[idx].imag() / (T) (kfft) * iiy / kw);
          psift2k[idx].imag(pField_k[idx].real() / (T) (kfft) * iiy / kw);
          psift3k[idx].real(-pField_k[idx].imag() / (T) (kfft) * iiz / kw);
          psift3k[idx].imag(pField_k[idx].real() / (T) (kfft) * iiz / kw);
        }
      }
    }

    psift1k[0] = complex<T>(0., 0.);
    psift2k[0] = complex<T>(0., 0.);
    psift3k[0] = complex<T>(0., 0.);

    fft(psift1k.data(), psift1k.data(), size,
        -1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
    fft(psift2k.data(), psift2k.data(), size, -1); //same
    fft(psift3k.data(), psift3k.data(), size, -1); //same

    pOff_x = std::make_shared<std::vector<T>>(size3, 0);
    pOff_y = std::make_shared<std::vector<T>>(size3, 0);
    pOff_z = std::make_shared<std::vector<T>>(size3, 0);

    //apply ZA:
#pragma omp parallel for schedule(static) default(shared) private(idx)
    for (size_t ix = 0; ix < size; ix++) {
      for (size_t iy = 0; iy < size; iy++) {
        for (size_t iz = 0; iz < size; iz++) {

          idx = (ix * size + iy) * size + iz;

          // position offset in physical coordinates
          (*pOff_x)[idx] = psift1k[idx].real();
          (*pOff_y)[idx] = psift2k[idx].real();
          (*pOff_z)[idx] = psift3k[idx].real();

        }
      }
    }

    cout << "done." << endl;

  }


  ///////////////////////////////////////
  //  Index manipulation routines
  ///////////////////////////////////////

  size_t getIndexFromIndexAndStep(size_t index, const Coordinate<int> &step) const {
    auto coord = getCellCoordinate(index);
    coord += step;
    return this->getCellIndex(coord); // N.B. does wrapping inside getIndex
  }

  void simWrap(T &x, T &y, T &z) const {
    x = fmod(x, simsize);
    if (x < 0) x += simsize;
    y = fmod(y, simsize);
    if (y < 0) y += simsize;
    z = fmod(z, simsize);
    if (z < 0) z += simsize;
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
  /*
  int xa=((int) floor((x0c-x0-dx/2)/dx));
  int ya=((int) floor((y0c-y0-dx/2)/dx));
  int za=((int) floor((z0c-z0-dx/2)/dx));
  return getIndexNoWrap(xa,ya,za);
  */

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
                    Coordinate<int>(xb, yb, zb)+1,
                    [&ids, this](const Coordinate<int> &cellCoord) {
                      ids.emplace_back(getCellIndex(cellCoord));
                    });

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
      pUnderlying->offsetLower.z, false),
    pUnderlying(pUnderlying) {

  }


  VirtualGrid(GridPtrType pUnderlying, T simsize, T gridsize,
              T dx, T x0, T y0, T z0, bool withField) :
    Grid<T>(simsize, gridsize, dx, x0, y0, z0, withField),
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

  bool pointsToGrid(Grid<T> *pOther) override {
    return pUnderlying.get() == pOther;
  }

  virtual bool fieldIsSuitableSize(const TRealField &field) override {
    return pUnderlying->fieldIsSuitableSize(field);
  }

  virtual bool fieldIsSuitableSize(const TField &field) override {
    return pUnderlying->fieldIsSuitableSize(field);
  }

  virtual T getFieldAt(size_t i, const TRealField &field) override {
    throw std::runtime_error("getFieldAt is not implemented for this type of VirtualGrid");
  }

  virtual complex<T> getFieldAt(size_t i, const TField &field) override {
    throw std::runtime_error("getFieldAt is not implemented for this type of VirtualGrid");
  }

  void gatherParticleList(std::vector<size_t> &targetArray) const override {
    pUnderlying->gatherParticleList(targetArray);
  }

  void distributeParticleList(const std::vector<size_t> &sourceArray) override {
    pUnderlying->distributeParticleList(sourceArray);
  }

  void clearParticleList() override {
    pUnderlying->clearParticleList();
  }

  size_t estimateParticleListSize() override {
    return pUnderlying->estimateParticleListSize();
  }

  virtual TField &getFieldFourier() override {
    return pUnderlying->getFieldFourier();
  }

  virtual TField &getFieldReal() override {
    return pUnderlying->getFieldReal();
  }

  virtual const TField &getField() const override {
    return pUnderlying->getField();
  }

  virtual bool isFieldFourier() const override {
    return pUnderlying->isFieldFourier();
  }

  virtual void zeldovich(T hfac, T particlecellMass) override {
    throw std::runtime_error("VirtualGrid - does not contain an actual field in memory");
  }


  virtual T getMass() const override {
    return pUnderlying->getMass();
  }

  void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {
    pUnderlying->getParticleNoWrap(id, x, y, z, vx, vy, vz, cellMassi, eps);
  }

  void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {
    pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
  }

  complex<T> getFieldAt(size_t i) override {
    return this->getFieldAt(i, this->pUnderlying->getFieldReal());
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
                   pUnderlying->offsetLower.z, false),
    factor(factor) {

    factor3 = factor * factor * factor;
  }

  virtual T getMass() const override {
    return this->pUnderlying->getMass() / factor3;
  }

  virtual void debugName(std::ostream &s) const override {
    s << "SuperSampleGrid";
  }

  void gatherParticleList(std::vector<size_t> &targetArray) const override {
    std::vector<size_t> underlyingArray;
    this->pUnderlying->gatherParticleList(underlyingArray);
    Grid<T>::upscaleParticleList(underlyingArray, targetArray, this->pUnderlying.get(), this);
  }

  void distributeParticleList(const std::vector<size_t> &sourceArray) override {
    std::vector<size_t> targetArray;
    Grid<T>::downscaleParticleList(sourceArray, targetArray, this, this->pUnderlying.get());
    this->pUnderlying->distributeParticleList(targetArray);
  }

  size_t estimateParticleListSize() override {
    return this->pUnderlying->estimateParticleListSize() * factor3;
  }

  virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    auto centroid = this->getCellCentroid(id);
    x = centroid.x;
    y = centroid.y;
    z = centroid.z;
    getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
  }

  void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {
    this->pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
    // adjust mass
    cellMassi /= factor3;
  }

  virtual T getFieldAt(size_t i, const TRealField &field) override {
    auto centroid = this->getCellCentroid(i);
    return this->pUnderlying->getFieldInterpolated(centroid, field);
  }

  virtual complex<T> getFieldAt(size_t i, const TField &field) override {
    auto centroid = this->getCellCentroid(i);
    return this->pUnderlying->getFieldInterpolated(centroid, field);
  }

};

template<typename T>
class OffsetGrid : public VirtualGrid<T> {
private:
  T xOffset, yOffset, zOffset;
  int xOffset_i, yOffset_i, zOffset_i;

protected:
  using typename Grid<T>::TField;
  using typename Grid<T>::TRealField;
  using typename Grid<T>::PtrTField;
  using typename Grid<T>::GridPtrType;


public:
  OffsetGrid(GridPtrType pUnderlying, T dx, T dy, T dz) :
    VirtualGrid<T>(pUnderlying),
    xOffset(dx), yOffset(dy), zOffset(dz) {

  }


  virtual void debugName(std::ostream &s) const override {
    s << "OffsetGrid";
  }

  virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    this->pUnderlying->getParticleNoWrap(id, x, y, z, vx, vy, vz, cellMassi, eps);
    x += xOffset;
    y += yOffset;
    z += zOffset;
    this->simWrap(x, y, z);
  }

  void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {
    x -= xOffset;
    y -= yOffset;
    z -= zOffset;
    this->simWrap(x, y, z);
    this->pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
    x += xOffset;
    y += yOffset;
    z += zOffset;
    this->simWrap(x, y, z);
  }

  void gatherParticleList(std::vector<size_t> &targetArray) const override {
    throw (std::runtime_error("gatherParticleList is not implemented for OffsetGrid"));
  }

  void distributeParticleList(const std::vector<size_t> &sourceArray) override {
    throw (std::runtime_error("distributeParticleList is not implemented for OffsetGrid"));
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
    int x, y, z;
    auto coord = this->getCellCoordinate(sec_id);
    coord += cellOffset;
    if (!this->pUnderlying->containsCell(coord))
      throw std::out_of_range("Out of range in SectionOfGrid::mapIndexToUnderlying");
    return this->pUnderlying->getCellIndex(coord);
  }

  size_t mapIndexFromUnderlying(size_t underlying_id) const {
    int x, y, z;
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
                   pUnderlying->offsetLower.z + deltaz * pUnderlying->dx, false),
    cellOffset(deltax, deltay, deltaz),
    posOffset(deltax * this->dx, deltay * this->dx, deltaz * this->dx),
    upperCell(cellOffset+pUnderlying->size)
  {

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

  virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    this->pUnderlying->getParticleNoWrap(mapIndexToUnderlying(id), x, y, z, vx, vy, vz, cellMassi, eps);
  }

  void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {

    this->pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);

  }

  virtual T getFieldAt(size_t i, const TRealField &field) override {
    return this->pUnderlying->getFieldAt(mapIndexToUnderlying(i), field);
  }

  virtual complex<T> getFieldAt(size_t i, const TField &field) override {
    return this->pUnderlying->getFieldAt(mapIndexToUnderlying(i), field);
  }

  void gatherParticleList(std::vector<size_t> &targetArray) const override {
    std::vector<size_t> underlyingArray;
    this->pUnderlying->gatherParticleList(underlyingArray);
    for (size_t ptcl: underlyingArray) {
      try {
        targetArray.push_back(this->mapIndexFromUnderlying(ptcl));
      } catch (std::out_of_range &e) {
        continue;
      }
    }
    std::sort(targetArray.begin(), targetArray.end());
  }

  void distributeParticleList(const std::vector<size_t> &sourceArray) override {
    throw (std::runtime_error("distributeParticleList is not implemented for OffsetGrid"));
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
                   pUnderlying->offsetLower.z, false),
    factor(factor) {
    //if(this->pUnderlying->size%factor!=0)
    //  throw std::runtime_error("SubSampleGrid - factor must be a divisor of the original grid size");

    factor3 = factor * factor * factor;
  }

  virtual void debugName(std::ostream &s) const override {
    s << "SubSampleGrid";
  }

  void gatherParticleList(std::vector<size_t> &targetArray) const override {
    std::vector<size_t> underlyingArray;
    this->pUnderlying->gatherParticleList(underlyingArray);
    Grid<T>::downscaleParticleList(underlyingArray, targetArray, this->pUnderlying.get(), this);
    // err << "SubSample gatherParticleList - underlying = " << underlyingArray.size() << " transformed = " <<targetArray.size() << endl;
  }

  void distributeParticleList(const std::vector<size_t> &sourceArray) override {
    std::vector<size_t> targetArray;
    Grid<T>::upscaleParticleList(sourceArray, targetArray, this, this->pUnderlying.get());
    this->pUnderlying->distributeParticleList(targetArray);
    // cerr << "SubSample distributeParticleList - source = " << sourceArray.size() << " transformed = " <<targetArray.size() << endl;
  }

  size_t estimateParticleListSize() override {
    return this->pUnderlying->estimateParticleListSize() / factor3;
  }


  virtual T getMass() const override {
    return this->pUnderlying->getMass() * factor3;
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

  virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    T xt, yt, zt, vxt, vyt, vzt, cellMassit, epst;

    x = 0;
    y = 0;
    z = 0;
    vx = 0;
    vy = 0;
    vz = 0;
    cellMassi = 0;
    eps = 0;

    int localFactor3 = forEachSubcell(id, [&](size_t sub_id) {
      this->pUnderlying->getParticleNoWrap(sub_id,
                                           xt, yt, zt, vxt, vyt, vzt, cellMassit, epst);
      x += xt;
      y += yt;
      z += zt;
      vx += vxt;
      vy += vyt;
      vz += vzt;
      eps += epst;
      cellMassi += cellMassit;
    });

    // most variables want an average, not a sum:
    x /= localFactor3;
    y /= localFactor3;
    z /= localFactor3;
    vx /= localFactor3;
    vy /= localFactor3;
    vz /= localFactor3;
    eps /= localFactor3;

    // cell mass wants to be a sum over the *entire* cell (even if
    // some subcells are missing):
    cellMassi *= factor3 / localFactor3;

  }

  virtual complex<T> getFieldAt(size_t i, const TField &field) {
    complex<T> returnVal(0);
    int localFactor3 = forEachSubcell(i, [this, &returnVal, &field](size_t local_id) {
      returnVal += field[local_id];
    });
    return returnVal / T(localFactor3);
  }

  virtual T getFieldAt(size_t i, const TRealField &field) {
    T returnVal(0);
    int localFactor3 = forEachSubcell(i, [this, &returnVal, &field](size_t local_id) {
      returnVal += field[local_id];
    });
    return returnVal / localFactor3;
  }


  void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {
    this->pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
    cellMassi *= factor3;
  }

};


template<typename T>
class MassScaledGrid : public VirtualGrid<T> {
protected:
  using typename Grid<T>::TField;
  using typename Grid<T>::TRealField;
  using typename Grid<T>::PtrTField;
  using typename Grid<T>::GridPtrType;

  T massFac;

public:
  MassScaledGrid(GridPtrType pUnderlying, T massFac) :
    VirtualGrid<T>(pUnderlying,
                   pUnderlying->simsize, pUnderlying->size,
                   pUnderlying->dx, pUnderlying->offsetLower.x, pUnderlying->offsetLower.y,
                   pUnderlying->offsetLower.z, false),
    massFac(massFac) {
  }

  virtual void debugName(std::ostream &s) const override {
    s << "MassScaledGrid";
  }

  virtual T getMass() const override {
    return this->pUnderlying->getMass() * massFac;
  }

  virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const {
    this->pUnderlying->getParticleNoWrap(id, x, y, z, vx, vy, vz, cellMassi, eps);
    cellMassi *= massFac;
  }

  void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override {
    this->pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
    cellMassi *= massFac;
  }

};


#endif

