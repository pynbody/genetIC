#ifndef IC_WINDOW_HPP
#define IC_WINDOW_HPP

#include "coordinate.hpp"
#include <cmath>

template<typename T>
class Window {
protected:
  T wrapLength;
  Coordinate<T> upperCornerExclusive;
  Coordinate<T> lowerCornerInclusive;


  bool withinWrapped(T lowInclusive, T highExclusive, T testVal) const {
    if (highExclusive >= lowInclusive)
      return testVal < highExclusive && testVal >= lowInclusive;
    else
      return testVal < highExclusive || testVal >= lowInclusive;
  }

  template<typename S>
  static S smallIncrement(S val) {
    return std::nextafter(val, std::numeric_limits<S>::infinity());
  }

  static int smallIncrement(int i) {
    return i + 1;
  }

  template<typename S>
  static S smallDecrement(S val) {
    return std::nextafter(val, -std::numeric_limits<S>::infinity());
  }

  static int smallDecrement(int i) {
    return i - 1;
  }

  static Coordinate<T> smallIncrement(const Coordinate<T> &input) {
    return Coordinate<T>(smallIncrement(input.x), smallIncrement(input.y), smallIncrement(input.z));
  }

  static Coordinate<T> smallDecrement(const Coordinate<T> &input) {
    return Coordinate<T>(smallDecrement(input.x), smallDecrement(input.y), smallDecrement(input.z));
  }


public:
  Window() : wrapLength(0) { }

  Window(T wrapLength, Coordinate<T> initialPosition) : wrapLength(wrapLength) {
    upperCornerExclusive = smallIncrement(initialPosition);
    lowerCornerInclusive = initialPosition;
    assert(this->contains(initialPosition));
  }

  Window(T wrapLength, Coordinate<T> lowerCornerInclusive, Coordinate<T> upperCornerExclusive) :
      wrapLength(wrapLength), upperCornerExclusive(upperCornerExclusive), lowerCornerInclusive(lowerCornerInclusive) {

  }

  Coordinate<T> getCurrentCentre() const {
    return wrap(lowerCornerInclusive + wrap(smallDecrement(upperCornerExclusive) - lowerCornerInclusive) / 2);
  }

  T wrap(T source) const {
    if (source < 0) source += wrapLength;
    if (source >= wrapLength) source -= wrapLength;
    return source;
  }

  //! Wrap a coordinate into the fundamental domain [0,wrapLength)
  Coordinate<T> wrap(const Coordinate<T> &source) const {
    Coordinate<T> result;
    result.x = wrap(source.x);
    result.y = wrap(source.y);
    result.z = wrap(source.z);
    return result;
  }


  //! Wrap a difference into the closest possible offset [-wrapLength/2, wrapLength/2)
  T getWrappedOffset(T a, T b) const {
    T distance = b - a;
    if (distance < -wrapLength / 2) distance += wrapLength;
    if (distance > wrapLength / 2) distance -= wrapLength;
    return distance;
  }

  Coordinate<T> getWrappedOffset(const Coordinate<T> &a, const Coordinate<T> &b) const {
    return Coordinate<T>(getWrappedOffset(a.x, b.x), getWrappedOffset(a.y, b.y), getWrappedOffset(a.z, b.z));
  }

  Coordinate<T> getSizes() const {
    Coordinate<T> size = upperCornerExclusive - lowerCornerInclusive;
    return wrap(size);
  }

  Coordinate<T> getLowerCornerInclusive() const {
    return lowerCornerInclusive;
  }

  T getMaximumDimension() const {
    auto offset = getSizes();
    return std::max({offset.x, offset.y, offset.z});
  }

  void expandToInclude(const Coordinate<T> &include, T Coordinate<T>::*coord) {
    if (!withinWrapped(lowerCornerInclusive.*coord, upperCornerExclusive.*coord, include.*coord)) {
      T sizeIfMovingLower = wrap(upperCornerExclusive.*coord - include.*coord - 1) + 1;
      T sizeIfMovingUpper = wrap(smallIncrement(include.*coord) - lowerCornerInclusive.*coord - 1) + 1;
      bool canMoveLower = withinWrapped(include.*coord, upperCornerExclusive.*coord, lowerCornerInclusive.*coord);
      bool canMoveUpper = withinWrapped(lowerCornerInclusive.*coord, smallIncrement(include.*coord),
                                        smallDecrement(upperCornerExclusive.*coord));

      assert(canMoveLower || canMoveUpper);

      if (canMoveLower && (sizeIfMovingLower < sizeIfMovingUpper || !canMoveUpper))
        lowerCornerInclusive.*coord = include.*coord;
      else
        upperCornerExclusive.*coord = smallIncrement(include.*coord);
    }
  }

  void expandToInclude(const Coordinate<T> &include) {
    expandToInclude(include, &Coordinate<T>::x);
    expandToInclude(include, &Coordinate<T>::y);
    expandToInclude(include, &Coordinate<T>::z);
  }

  void expandSymmetricallyToSize(T newSize) {
    assert(getMaximumDimension() <= newSize);
    auto centre = getCurrentCentre();

    auto oldLCI = lowerCornerInclusive;
    auto oldUCE = upperCornerExclusive;

    lowerCornerInclusive = wrap(centre - newSize / 2);
    if (newSize % 2 == 0)
      upperCornerExclusive = wrap(centre + newSize / 2);
    else
      upperCornerExclusive = wrap(centre + (newSize / 2 + 1));

    // There is an edge case where T is int, newSize is odd, and also getSize()==newSize, where the new
    // window might be aligned too far up to capture the old lower corner. Fix this:
    if (lowerCornerInclusive.x > oldLCI.x) {
      lowerCornerInclusive.x -= 1;
      upperCornerExclusive.x -= 1;
    }
    if (lowerCornerInclusive.y > oldLCI.z) {
      lowerCornerInclusive.y -= 1;
      upperCornerExclusive.y -= 1;
    }
    if (lowerCornerInclusive.z > oldLCI.z) {
      lowerCornerInclusive.z -= 1;
      upperCornerExclusive.z -= 1;
    }

    assert(this->contains(oldLCI));
    assert(this->contains(oldUCE - 1));
    assert(getSizes() == newSize);

  }

  bool contains(const Coordinate<T> &test) const  {
    bool inX, inY, inZ;
    inX = withinWrapped(lowerCornerInclusive.x, upperCornerExclusive.x, test.x);
    inY = withinWrapped(lowerCornerInclusive.y, upperCornerExclusive.y, test.y);
    inZ = withinWrapped(lowerCornerInclusive.z, upperCornerExclusive.z, test.z);
    return inX && inY && inZ;
  }

  bool containsWithBorderSafety(const Coordinate<T> &test, T safety) const {
    bool inX, inY, inZ;
    inX = withinWrapped(lowerCornerInclusive.x + safety, upperCornerExclusive.x - safety, test.x);
    inY = withinWrapped(lowerCornerInclusive.y + safety, upperCornerExclusive.y - safety, test.y);
    inZ = withinWrapped(lowerCornerInclusive.z + safety, upperCornerExclusive.z - safety, test.z);
    return inX && inY && inZ;
  }

};

#endif //IC_WINDOW_HPP
