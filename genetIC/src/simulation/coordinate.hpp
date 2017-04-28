#ifndef IC_COORDINATE_HPP
#define IC_COORDINATE_HPP

#include <cmath>
#include <functional>

/*!
    \class Coordinate
    \brief Definition of coordinates for grid cells and particles.

 */
template<typename T>
class Coordinate {
public:
  T x, y, z;

  Coordinate(T x, T y, T z) : x(x), y(y), z(z) {
  }

  Coordinate(T *vec) : x(vec[0]), y(vec[1]), z(vec[2]) {
  }

  Coordinate(T value) : x(value), y(value), z(value) {}

  Coordinate() {

  }

  template<typename S>
  Coordinate(const Coordinate<S> &other) : x(other.x), y(other.y), z(other.z) {}


  Coordinate<T> operator+(const Coordinate<T> &other) const {
    return Coordinate<T>(x + other.x, y + other.y, z + other.z);
  }

  Coordinate<T> operator+(T offset) const {
    return Coordinate<T>(x + offset, y + offset, z + offset);
  }

  Coordinate<T> operator-(const Coordinate<T> &other) const {
    return Coordinate<T>(x - other.x, y - other.y, z - other.z);
  }

  Coordinate<T> operator-(T offset) const {
    return Coordinate<T>(x - offset, y - offset, z - offset);
  }

  Coordinate<T> operator/(T factor) const {
    return Coordinate<T>(x / factor, y / factor, z / factor);
  }

  Coordinate<T> operator*(T factor) const {
    return Coordinate<T>(x * factor, y * factor, z * factor);
  }

  operator std::tuple<T &, T &, T &>() const {
    //return std::make_tuple(const_cast<T&>(x),const_cast<T&>(y),const_cast<T&>(z));
    return std::tuple<T &, T &, T &>(const_cast<T &>(x),
                                     const_cast<T &>(y),
                                     const_cast<T &>(z));
  };

  void operator+=(const Coordinate<T> &other) {
    x += other.x;
    y += other.y;
    z += other.z;
  }

  void operator+=(T offset) {
    x += offset;
    y += offset;
    z += offset;
  }

  void operator-=(const Coordinate<T> &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
  }

  void operator-=(T offset) {
    x -= offset;
    y -= offset;
    z -= offset;
  }

  void operator/=(T factor) {
    x /= factor;
    y /= factor;
    z /= factor;
  }

  void operator*=(T factor) {
    x *= factor;
    y *= factor;
    z *= factor;
  }

  bool operator==(const Coordinate<T> &other) const {
    return x == other.x && y == other.y && z == other.z;
  }

  bool operator!=(const Coordinate<T> &other) const {
    return !((*this) == other);
  }

  bool almostEqual(const Coordinate<T> &other, T tol=1e-10) const {
    return abs(x-other.x)<tol && abs(y-other.y)<tol && abs(z-other.z)<tol;
  }


};

template<typename T>
void iterateOverCube(const Coordinate<T> &lowerCornerInclusive, const Coordinate<T> &upperCornerExclusive,
                     std::function<void(const Coordinate<T> &)> callback) {
  T x, y, z;
  for (x = lowerCornerInclusive.x; x < upperCornerExclusive.x; ++x) {
    for (y = lowerCornerInclusive.y; y < upperCornerExclusive.y; ++y) {
      for (z = lowerCornerInclusive.z; z < upperCornerExclusive.z; ++z) {
        callback(Coordinate<T>(x, y, z));
      }
    }
  }
}

template<typename T>
Coordinate<int> floor(const Coordinate<T> &coord) {
  return Coordinate<int>(int(std::floor(coord.x)),
                         int(std::floor(coord.y)),
                         int(std::floor(coord.z)));

}

#endif //IC_COORDINATE_HPP
