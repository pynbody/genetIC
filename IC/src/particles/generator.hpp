//
// Created by Andrew Pontzen on 20/12/2016.
//

#ifndef IC_GENERATOR_HPP
#define IC_GENERATOR_HPP


#include <complex>
#include "src/field.hpp"
#include "src/cosmo.hpp"

namespace particle {

  template<typename T>
  class ParticleGenerator : public std::enable_shared_from_this<ParticleGenerator<T>> {
  protected:
    const Grid <T> &grid;
    ParticleGenerator(const Grid<T> &grid) : grid(grid) {

    }

  public:

    particle::Particle<T> getParticle(size_t id) const {
      particle::Particle<T> particle = getParticleNoWrap(id);
      grid.simWrap(particle.pos);
      return particle;
    }

    const Grid<T> & getGrid() const {
      return grid;
    }

    virtual void recalculate() = 0;
    virtual particle::Particle<T> getParticleNoWrap(size_t id) const =0;
    virtual T getMass() const =0;
    virtual T getEps() const =0;
    virtual particle::Particle<T> getParticleNoOffset(size_t id) const =0;
    virtual void getParticleFromOffset(particle::Particle<T> &particle) const =0;


  };


}
#endif //IC_GENERATOR_HPP
