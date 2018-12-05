#ifndef IC_PARTICLE_HPP
#define IC_PARTICLE_HPP

#include "src/simulation/coordinate.hpp"

/*!
    \namespace particle
    \brief Generate particles for N-body zoomed simulations from multiple levels grids.
*/
namespace particle {
  /*! \class Particle
    \brief Defines a particle, and stores necessary properties (position, velocity, mass, and cell softening scale)
  */
  template<typename T>
  class Particle {
  public:
    Coordinate<T> pos; //!< Position in comoving co-ordinates
    Coordinate<T> vel; //!< Velcoity in comoving co-ordinates
    T mass; //!< Mass of the particle
    T soft; //!< Cell softening scale of this particle

    //! Constructor which initialises everything to 0
    Particle(int) : pos(T(0)), vel(T(0)), mass(0), soft(0) {

    }

    //! Constructor leaving everything undefined.
    Particle() {}

  };
}

#endif
