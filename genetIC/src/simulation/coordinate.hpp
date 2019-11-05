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
  T x, y, z; //!< x,y,and z positions defining the co-ordinate

  //! Creates a co-ordinate from specified floats
  Coordinate(T x, T y, T z) : x(x), y(y), z(z) {
  }

  //! Creates a co-ordinate from an array of floats (doesn't check bounds!)
  Coordinate(T *vec) : x(vec[0]), y(vec[1]), z(vec[2]) {
  }

  //! Creates a co-ordinate with all three components the same
  Coordinate(T value) : x(value), y(value), z(value) {}

  //! Creates a co-ordinate with zero in all components
  Coordinate() : x(0), y(0), z(0) {}

  //! Creates a co-ordinate with array initialisation
  template<size_t N>
  Coordinate(const T(&vals)[N]) {
    static_assert(N == 3, "Coordinate can only be initialized with three elements");
    x = vals[0];
    y = vals[1];
    z = vals[2];
  }

  //! Creates a co-ordinate from another co-ordinates
  template<typename S>
  explicit Coordinate(const Coordinate<S> &other) : x(other.x), y(other.y), z(other.z) {}


  //! Returns another co-ordinate added to this one
  Coordinate<T> operator+(const Coordinate<T> &other) const {
    return Coordinate<T>(x + other.x, y + other.y, z + other.z);
  }

  //! Returns a co-ordinate with a constant added to all components of this one
  Coordinate<T> operator+(T offset) const {
    return Coordinate<T>(x + offset, y + offset, z + offset);
  }

  //! Subtracts a co-ordinate from this one
  Coordinate<T> operator-(const Coordinate<T> &other) const {
    return Coordinate<T>(x - other.x, y - other.y, z - other.z);
  }

  //! Subtracts a constant from all components of this one
  Coordinate<T> operator-(T offset) const {
    return Coordinate<T>(x - offset, y - offset, z - offset);
  }

  //! Returns the co-ordinate with all components flipped by a sign
  Coordinate<T> operator-() const {
    return Coordinate<T>(-x, -y, -z);
  }

  //! Returns the co-ordinate divided by some factor
  Coordinate<T> operator/(T factor) const {
    return Coordinate<T>(x / factor, y / factor, z / factor);
  }

  //! Returns the co-ordinate mutliplied by some factor
  Coordinate<T> operator*(T factor) const {
    return Coordinate<T>(x * factor, y * factor, z * factor);
  }

  //! Converts the co-ordinate to sts::tuple
  operator std::tuple<T &, T &, T &>() const {
    return std::tuple<T &, T &, T &>(const_cast<T &>(x),
                                     const_cast<T &>(y),
                                     const_cast<T &>(z));
  };

  //! Adds some co-ordinate to this one.
  void operator+=(const Coordinate<T> &other) {
    x += other.x;
    y += other.y;
    z += other.z;
  }

  //! Adds a constant offset to all components of this co-ordinate
  void operator+=(T offset) {
    x += offset;
    y += offset;
    z += offset;
  }

  //! Subtracts another co-ordinate from all this co-ordinate
  void operator-=(const Coordinate<T> &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
  }

  //! Subtracts a constant offset from all components for this co-ordinate
  void operator-=(T offset) {
    x -= offset;
    y -= offset;
    z -= offset;
  }

  //! Divides this co-ordinate by the specified factor
  void operator/=(T factor) {
    x /= factor;
    y /= factor;
    z /= factor;
  }

  //! Multiplies this co-ordinate by the specified factor
  void operator*=(T factor) {
    x *= factor;
    y *= factor;
    z *= factor;
  }

  T & operator[](size_t i) {
    if(i==0)
      return x;
    else if(i==1)
      return y;
    else if(i==2)
      return z;
    else
      throw std::out_of_range("Coordinates always have three dimensions");
  }

  //! Returns true if all components of the other co-ordinate match this one
  bool operator==(const Coordinate<T> &other) const {
    return x == other.x && y == other.y && z == other.z;
  }

  //! Returns true if at least one component of the other co-ordinate differs from this one
  bool operator!=(const Coordinate<T> &other) const {
    return !((*this) == other);
  }

  //! Returns true if the other co-ordinate is 'almost' equal to this one, within the specified tolerance
  bool almostEqual(const Coordinate<T> &other, T tol = 1e-10) const {
    return std::abs(x - other.x) < tol && std::abs(y - other.y) < tol && std::abs(z - other.z) < tol;
  }

  /*! \brief Returns true if the co-ordinate points to somewhere that lies in the window  defined by the two specified co-ordinates
    \param lowerCornerInclusive - lower left front corner of the window, including that point
    \param upperCornerExclusive - upper right back corner of the window, excluding that point
  */
  bool inWindow(const Coordinate<T> &lowerCornerInclusive, const Coordinate<T> &upperCornerExclusive) const {
    return x >= lowerCornerInclusive.x && x < upperCornerExclusive.x &&
           y >= lowerCornerInclusive.y && y < upperCornerExclusive.y &&
           z >= lowerCornerInclusive.z && z < upperCornerExclusive.z;
  }


};

/*! \brief Iterates over all points in the window defined by the two co-ordinates, calling a function at each point.
    \param lowerCornerInclusive - lower left front corner of the window, including that point
    \param upperCornerExclusive - upper right back corner of the window, excluding that point
    \param callback - function to apply to each point in the window
*/
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

//! Floors the components of an integer co-ordinate
template<typename T>
Coordinate<int> floor(const Coordinate<T> &coord) {
  return Coordinate<int>(int(std::floor(coord.x)),
                         int(std::floor(coord.y)),
                         int(std::floor(coord.z)));

}

//! Rounds the components of a specified co-ordinate to the nearest integer, and converts the type from T to OutputType
template<typename OutputType, typename T>
Coordinate<OutputType> round(const Coordinate<T> &coord) {
  return Coordinate<OutputType>(OutputType(std::round(coord.x)),
                                OutputType(std::round(coord.y)),
                                OutputType(std::round(coord.z)));

}

//! Outputs the components of a co-ordinate to the specified stream
template<typename T>
std::ostream &operator<<(std::ostream &stream, const Coordinate<T> &coord) {
  stream << "(" << coord.x << ", " << coord.y << ", " << coord.z << ")";
  return stream;
}

#endif //IC_COORDINATE_HPP
