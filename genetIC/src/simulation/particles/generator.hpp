#ifndef IC_GENERATOR_HPP
#define IC_GENERATOR_HPP


#include <complex>
#include "src/simulation/field/field.hpp"
#include "src/simulation/particles/particle.hpp"

namespace particle {

  /** A class to evaluate the results of the particle generation on a given grid. This allows for
   * flexible interpolation etc when the grid being evaluated on is not necessarily the same as the particles
   * were actually generated on.
   */
  template<typename GT>
  class ParticleEvaluator : public std::enable_shared_from_this<ParticleEvaluator<GT>> {
  protected:
    using T = tools::datatypes::strip_complex<GT>;
    const grids::Grid<T> &grid;

    ParticleEvaluator(const grids::Grid<T> &grid) : grid(grid) { };

  public:

    virtual Particle <T> getParticleNoWrap(size_t id) const =0;

    virtual T getMass() const =0;

    virtual T getEps() const =0;

    virtual Particle <T> getParticleNoOffset(size_t id) const =0;

    const grids::Grid<T> &getGrid() const {
      return grid;
    }

    Particle <T> getParticle(size_t id) const {
      Particle<T> particle = getParticleNoWrap(id);
      grid.simWrap(particle.pos);
      return particle;
    }
  };


  /** A class to generate particles starting from cells of grids. Note this is a pure virtual class because
   * we require a specific implementation (currently Zeldovich is the only option!)
   */
  template<typename GT>
  class ParticleGenerator : public std::enable_shared_from_this<ParticleGenerator<GT>> {
  public:
    using T = tools::datatypes::strip_complex<GT>;

  protected:
    const grids::Grid<T> &grid;

    ParticleGenerator(const grids::Grid<T> &grid) : grid(grid) {

    }

  public:
    const grids::Grid<T> &getGrid() const {
      return grid;
    }

    virtual void recalculate() = 0;

    virtual std::shared_ptr<ParticleEvaluator<GT>> getEvaluator(const grids::Grid<T> &forGrid) const =0;

  };




}
#endif
