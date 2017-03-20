#ifndef IC_PARTICLE_HPP
#define IC_PARTICLE_HPP

#include "src/simulation/coordinate.hpp"

/*!
    \namespace particle
    \brief Generate particles for N-body zoomed simulations from multiple levels grids.
*/
namespace particle {
  template<typename T>
  class Particle {
  public:
    Coordinate<T> pos;
    Coordinate<T> vel;
    T mass;
    T soft;

    Particle(int) : pos(T(0)), vel(T(0)), mass(0), soft(0) {

    }

    Particle() {}

  };
}

#endif //IC_PARTICLE_HPP
