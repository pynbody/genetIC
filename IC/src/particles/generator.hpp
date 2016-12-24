//
// Created by Andrew Pontzen on 20/12/2016.
//

#ifndef IC_GENERATOR_HPP
#define IC_GENERATOR_HPP


#include <complex>
#include "src/field.hpp"
#include "src/cosmo.hpp"
#include "particle.hpp"

namespace particle {

  template<typename T>
  class Particle;

  template<typename T>
  class ParticleGenerator : public std::enable_shared_from_this<ParticleGenerator<T>> {
  protected:
    const Grid <T> &grid;
    ParticleGenerator(const Grid<T> &grid) : grid(grid) {

    }

  public:

    Particle<T> getParticle(const Grid<T> & onGrid, size_t id) const {
      Particle<T> particle = getParticleNoWrap(onGrid, id);
      grid.simWrap(particle.pos);
      return particle;
    }

    const Grid<T> & getGrid() const {
      return grid;
    }

    virtual void recalculate() = 0;
    virtual Particle<T> getParticleNoWrap(const Grid<T> & onGrid, size_t id) const =0;
    virtual T getMass(const Grid<T> & onGrid) const =0;
    virtual T getEps(const Grid<T> & onGrid) const =0;
    virtual Particle<T> getParticleNoOffset(const Grid<T> & onGrid, size_t id) const =0;
    // virtual void getParticleFromOffset(const Grid<T> & onGrid, particle::Particle<T> &particle) const =0;


  };


}
#endif //IC_GENERATOR_HPP
