//
// Created by Andrew Pontzen on 20/12/2016.
//

#ifndef IC_GENERATOR_HPP
#define IC_GENERATOR_HPP


#include <complex>
#include "src/simulation/field/field.hpp"
#include "src/cosmology/parameters.hpp"
#include "src/simulation/particles/particle.hpp"

namespace particle {

  template<typename T>
  class Particle;

  /** A class to generate particles starting from cells of grids. Note this is a pure virtual class because
   * we require a specific implementation (currently Zeldovich is the only option!)
   */
  template<typename GT>
  class ParticleGenerator : public std::enable_shared_from_this<ParticleGenerator<GT>> {
  public:
    using T = tools::datatypes::strip_complex<GT>;

  protected:
    const grids::Grid <T> &grid;
    ParticleGenerator(const grids::Grid<T> &grid) : grid(grid) {

    }

  public:

    Particle<T> getParticle(const grids::Grid<T> & onGrid, size_t id) const {
      Particle<T> particle = getParticleNoWrap(onGrid, id);
      grid.simWrap(particle.pos);
      return particle;
    }

    const grids::Grid<T> & getGrid() const {
      return grid;
    }

    virtual void recalculate() = 0;
    virtual Particle<T> getParticleNoWrap(const grids::Grid<T> & onGrid, size_t id) const =0;
    virtual T getMass(const grids::Grid<T> & onGrid) const =0;
    virtual T getEps(const grids::Grid<T> & onGrid) const =0;
    virtual Particle<T> getParticleNoOffset(const grids::Grid<T> & onGrid, size_t id) const =0;
    // virtual void getParticleFromOffset(const grids::Grid<T> & onGrid, particle::Particle<T> &particle) const =0;


  };


}
#endif //IC_GENERATOR_HPP
