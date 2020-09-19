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


  //! Wrap a difference into the closest possible offset [-wrapLength/2, wrapLength/2)
  T getWrappedOffset(T a, T b) const {
    T distance = b - a;
    if (distance < -wrapLength / 2) distance += wrapLength;
    if (distance > wrapLength / 2) distance -= wrapLength;
    return distance;
  }

  //! Get the offset between two co-ordinates, taking into account wrapping
  Coordinate<T> getWrappedOffset(const Coordinate<T> &a, const Coordinate<T> &b) const {
    return Coordinate<T>(getWrappedOffset(a.x, b.x), getWrappedOffset(a.y, b.y), getWrappedOffset(a.z, b.z));
  }

  //! Returns a co-ordinate pointing from the lower front left to the upper back right hand corner, wrapped into the fundamental domain
  Coordinate<T> getSizes() const {
    Coordinate<T> size = upperCornerExclusive - lowerCornerInclusive;
    return wrap(size);
  }

  //! Returns the position of the lower front left hand corner
  Coordinate<T> getLowerCornerInclusive() const {
    return lowerCornerInclusive;
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
    auto centre = getCurrentCentre();

    auto oldLCI = lowerCornerInclusive;
    auto oldUCE = upperCornerExclusive;

    lowerCornerInclusive = wrap(centre - newSize / 2);
    if (newSize % 2 == 0)
      upperCornerExclusive = wrap(centre + newSize / 2 - 1) + 1;
    else
      upperCornerExclusive = wrap(centre + newSize / 2 ) + 1;

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

    assert(this->contains(oldLCI));
    assert(this->contains(oldUCE - 1));

    assert(getSizes() == wrap(newSize));

  }

  //! Moves the window into the fundamental wrapping domain
  void clampToFundamentalDomain() {
    auto oldSize = wrap(upperCornerExclusive - lowerCornerInclusive);
    if(oldSize.x>=wrapLength || oldSize.y>=wrapLength || oldSize.z>=wrapLength) {
      throw std::runtime_error("The window is too large to be moved into the fundamental domain");
      // Hopefully a user never encounters the above error. If this error is reached, the calling code should
      // be examined for inconsistencies -- why is it trying to shift something into a box it doesn't fit in?
    }

    if (lowerCornerInclusive.x<=0) {
      upperCornerExclusive.x-=lowerCornerInclusive.x;
      lowerCornerInclusive.x=0;
    }
    if (lowerCornerInclusive.y<=0) {
      upperCornerExclusive.y-=lowerCornerInclusive.y;
      lowerCornerInclusive.y=0;
    }
    if (lowerCornerInclusive.z<=0) {
      upperCornerExclusive.z-=lowerCornerInclusive.z;
      lowerCornerInclusive.z=0;
    }

    if (upperCornerExclusive.x>wrapLength) {
      lowerCornerInclusive.x = wrapLength - oldSize.x;
      upperCornerExclusive.x=wrapLength;
    }
    if (upperCornerExclusive.y>wrapLength) {
      lowerCornerInclusive.y = wrapLength - oldSize.y;
      upperCornerExclusive.y=wrapLength;
    }
    if (upperCornerExclusive.z>wrapLength) {
      lowerCornerInclusive.z = wrapLength - oldSize.z;
      upperCornerExclusive.z=wrapLength;
    }

    auto newSize = upperCornerExclusive - lowerCornerInclusive;
    assert(oldSize==newSize);
  }

  //! Returns true if the window contains the test co-ordinate
  bool contains(const Coordinate<T> &test) const {
    bool inX, inY, inZ;
    inX = withinWrapped(lowerCornerInclusive.x, upperCornerExclusive.x, test.x);
    inY = withinWrapped(lowerCornerInclusive.y, upperCornerExclusive.y, test.y);
    inZ = withinWrapped(lowerCornerInclusive.z, upperCornerExclusive.z, test.z);
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
