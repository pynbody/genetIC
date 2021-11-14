#ifndef IC_WINDOW_HPP
#define IC_WINDOW_HPP

#include "coordinate.hpp"
#include <cmath>

/*! \class Window
    \brief Defines a window, ie, a cuboid region defined by a lower front left hand side co-ordinate, and an upper back right hand side co-ordinate
*/
template<typename T>
class Window {
protected:
  T wrapLength; //!< Periodicity scale after which co-ordinates wrap back on themselves
  Coordinate<T> upperCornerExclusive; //!< Upper back right hand side corner, but not included in the window itself
  Coordinate<T> lowerCornerInclusive; //!< Lower front left hand side corner, included in the window


  /*! \brief Returns true if a value is in a given range, taking into account periodicity
      \param lowInclusive - lower bound of the range (included in range)
      \param highExclusive - upper bound of the range (excluded from range)
      \param testVal - value we wish to check
  */
  bool withinWrapped(T lowInclusive, T highExclusive, T testVal) const {
    if (highExclusive >= lowInclusive)
      return testVal < highExclusive && testVal >= lowInclusive;
    else
      return testVal < highExclusive || testVal >= lowInclusive;
  }

  //! Increments a variable by the smallest amount computationally possible
  template<typename S>
  static S smallIncrement(S val) {
    return std::nextafter(val, std::numeric_limits<S>::infinity());
  }

  //! Specialisation to integers
  static int smallIncrement(int i) {
    return i + 1;
  }

  //! Decrements a variable by the smallest amount computationally possible
  template<typename S>
  static S smallDecrement(S val) {
    return std::nextafter(val, -std::numeric_limits<S>::infinity());
  }

  //! Specialisation to integers
  static int smallDecrement(int i) {
    return i - 1;
  }

  //! Performs a small increment for all components of a co-ordinate
  static Coordinate<T> smallIncrement(const Coordinate<T> &input) {
    return Coordinate<T>(smallIncrement(input.x), smallIncrement(input.y), smallIncrement(input.z));
  }

  //! Performs a small decrement for all components of a co-ordinate
  static Coordinate<T> smallDecrement(const Coordinate<T> &input) {
    return Coordinate<T>(smallDecrement(input.x), smallDecrement(input.y), smallDecrement(input.z));
  }


public:
  //! Construct an empty window, with zero wrap length
  Window() : wrapLength(0) {}

  /*! \brief Constructs a window with the given wrap length, and the smallest possible size around the specified initial position
      \param wrapLength - periodicity scale
      \param initialPosition - lower front left hand corner of window
  */
  Window(T wrapLength, Coordinate<T> initialPosition) : wrapLength(wrapLength) {
    upperCornerExclusive = smallIncrement(initialPosition);
    lowerCornerInclusive = initialPosition;
    assert(this->contains(initialPosition));
  }

  /*! \brief Constructs a window with the specified wrap length, and lower/upper corner
      \param wrapLength - periodicity scale
      \param lowInclusive - lower bound of the range (included in range)
      \param highExclusive - upper bound of the range (excluded from range)
  */
  Window(T wrapLength, Coordinate<T> lowerCornerInclusive, Coordinate<T> upperCornerExclusive) :
    wrapLength(wrapLength), upperCornerExclusive(upperCornerExclusive), lowerCornerInclusive(lowerCornerInclusive) {
    wrapCorners();
  }

  //! Returns a co-ordinate pointing to the centre of the current window
  Coordinate<T> getCurrentCentre() const {
    return wrap(lowerCornerInclusive + wrap(smallDecrement(upperCornerExclusive) - lowerCornerInclusive) / 2);
  }

  //! Wraps the specified variable with the defined periodicity, so it lies in the fundamental domain
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

  void wrapCorners() {
    this->lowerCornerInclusive = wrap(lowerCornerInclusive);

    // because the upper corner is exclusive, the actual point to wrap is the one just inside
    // This is an important distinction if the window covers the whole domain (since otherwise
    // the upper corner can wrap back round to zero, making a full-size window look like a
    // zero-size window!)
    this->upperCornerExclusive = smallIncrement(wrap(smallDecrement(upperCornerExclusive)));
  }

  void operator-=(const Coordinate<T> &shift) {
    lowerCornerInclusive-=shift;
    upperCornerExclusive-=shift;
    wrapCorners();
  }

  void operator+=(const Coordinate<T> &shift) {
    lowerCornerInclusive+=shift;
    upperCornerExclusive+=shift;
    wrapCorners();
  }


  //! Wrap a difference into the closest possible offset [-wrapLength/2, wrapLength/2)
  T getWrappedOffset(T a, T b) const {
    T distance = b - a;
    if (distance < -wrapLength / 2) distance += wrapLength;
    if (distance > wrapLength / 2) distance -= wrapLength;
    return distance;
  }

  //! Get the offset between two co-ordinates, b-a, taking into account wrapping
  Coordinate<T> getWrappedOffset(const Coordinate<T> &a, const Coordinate<T> &b) const {
    return Coordinate<T>(getWrappedOffset(a.x, b.x), getWrappedOffset(a.y, b.y), getWrappedOffset(a.z, b.z));
  }

  //! Returns a co-ordinate pointing from the lower front left to the upper back right hand corner, wrapped into the fundamental domain
  Coordinate<T> getSizes() const {
    Coordinate<T> size = upperCornerExclusive - lowerCornerInclusive;
    return smallIncrement(wrap(smallDecrement(size)));
  }

  //! Returns the position of the lower front left hand corner
  Coordinate<T> getLowerCornerInclusive() const {
    return lowerCornerInclusive;
  }

  //! Returns the integer cell window corresponding to this
  Window<int> convertPointToCell(T cellSize) const {
    Coordinate<int> windowLCI;
    Coordinate<int> windowUCE;
    for(int dir=0; dir<3; dir++) {
      windowLCI[dir] = tools::getRatioAndAssertPositiveInteger(lowerCornerInclusive[dir], cellSize);
      windowUCE[dir] = tools::getRatioAndAssertPositiveInteger(upperCornerExclusive[dir], cellSize);
    }
    int windowWrapLength = tools::getRatioAndAssertPositiveInteger(this->wrapLength, cellSize);
    return {windowWrapLength, windowLCI, windowUCE};
  }

  //! Returns the largest dimension of the cuboid window
  T getMaximumDimension() const {
    auto offset = getSizes();
    return std::max({offset.x, offset.y, offset.z});
  }

  /*! \brief Expand the window to include the specidied co-ordinate, in the specified direction only
    \param include - co-ordinate we need to include
    \param coord - x,y, or z, ie, direction to expand in.
  */
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

  //! Expand along all dimensions to include the specified co-ordinate in the window
  void expandToInclude(const Coordinate<T> &include) {
    expandToInclude(include, &Coordinate<T>::x);
    expandToInclude(include, &Coordinate<T>::y);
    expandToInclude(include, &Coordinate<T>::z);
  }

  //! Expand symmetrically in all directions to make the new distance between upper and lower corners equal to newSize
  void expandSymmetricallyToSize(T newSize) {
    assert(getMaximumDimension() <= newSize);
    auto oldLCI = lowerCornerInclusive;
    auto oldUCE = upperCornerExclusive;

    if(newSize==wrapLength) {
      lowerCornerInclusive = {0,0,0};
      upperCornerExclusive = {newSize, newSize, newSize};
    } else {
      auto centre = getCurrentCentre();


      lowerCornerInclusive = wrap(centre - newSize / 2);
      if (newSize % 2 == 0)
        upperCornerExclusive = wrap(centre + newSize / 2 - 1) + 1;
      else
        upperCornerExclusive = wrap(centre + newSize / 2) + 1;

      // Note in the above, because the corner is *exclusive*, the +1 happens after the wrapping

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
    }

    assert(this->contains(oldLCI));
    assert(this->contains(oldUCE - 1));

    assert(getSizes() == smallIncrement(wrap(smallDecrement(newSize))));

  }

  //! Move this window minimally so that it fits inside the other
  void moveInto(const Window<T> &other) {
    auto otherSizes = other.getSizes();
    auto thisSizes = getSizes();
    assert(otherSizes.x>=thisSizes.x && otherSizes.y>=thisSizes.y && otherSizes.z>=thisSizes.z);
    // If above assertion fails, code is trying to move this window into another window that it doesn't
    // actually fit inside! Look for bugs in the caller.

    assert(other.wrapLength == this->wrapLength);
    // If this fails, the wrapping length is ambiguous. The caller is asking for something dodgy.

    Coordinate<T> moveDistance;

    for(int dim=0; dim<3; dim++) {
      if (!other.contains(this->lowerCornerInclusive[dim], dim))
        moveDistance[dim] = getWrappedOffset(this->lowerCornerInclusive[dim], other.lowerCornerInclusive[dim]);

      if (!other.contains(smallDecrement(this->upperCornerExclusive[dim]), dim)) {
        assert(moveDistance[dim]==0); // if this fails, we are trying to move left already... we can't also move right!
        // suggests the window is too big (which should have been tested above already, so we shouldn't reach this point.)
        moveDistance[dim] = getWrappedOffset(this->upperCornerExclusive[dim], other.upperCornerExclusive[dim]);
      }

    }

    (*this)+=moveDistance;

    assert(other.contains(lowerCornerInclusive));
    assert(other.contains(smallDecrement(upperCornerExclusive)));


  }

  //! Returns true if the window contains the test co-ordinate in the specified dimension
  bool contains(T coord, int dim) const {
    return withinWrapped(lowerCornerInclusive[dim], upperCornerExclusive[dim], coord);
  }

  //! Returns true if the window contains the test co-ordinate
  bool contains(const Coordinate<T> &test) const {
    bool inX, inY, inZ;
    inX = contains(test.x, 0);
    inY = contains(test.y, 1);
    inZ = contains(test.z, 2);
    return inX && inY && inZ;
  }

  /*! \brief Returns true if the test co-ordinate lies within the window, and is not within safety threshold of the edge for any component
      \param test - co-ordinate to check
      \param safety - threshold - fails test if some component is closer to the edge than this
  */
  bool containsWithBorderSafety(const Coordinate<T> &test, T safety) const {
    bool inX, inY, inZ;
    inX = withinWrapped(lowerCornerInclusive.x + safety, upperCornerExclusive.x - safety, test.x);
    inY = withinWrapped(lowerCornerInclusive.y + safety, upperCornerExclusive.y - safety, test.y);
    inZ = withinWrapped(lowerCornerInclusive.z + safety, upperCornerExclusive.z - safety, test.z);
    return inX && inY && inZ;
  }

};

#endif //IC_WINDOW_HPP
