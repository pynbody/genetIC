#ifndef IC_GENERATOR_HPP
#define IC_GENERATOR_HPP


#include <complex>
#include "src/simulation/field/field.hpp"
#include "src/simulation/particles/particle.hpp"

namespace particle {

  /*! \class ParticleEvaluator
      \brief A class to evaluate the results of the particle generation on a given grid.

     This allows for flexible interpolation etc when the grid being evaluated on is not
     necessarily the same as the particles were actually generated on.
   */
  template<typename GT>
  class ParticleEvaluator : public std::enable_shared_from_this<ParticleEvaluator<GT>> {
  protected:
    using T = tools::datatypes::strip_complex<GT>;
    const grids::Grid<T> &grid; //!< Grid on which we are evaluating particles

    //! Constructor, which requires only the grid we want to evaluate on
    ParticleEvaluator(const grids::Grid<T> &grid) : grid(grid) {};

  public:

    //! Returns the particle at index id on the grid, without wrapping
    virtual Particle <T> getParticleNoWrap(size_t id) const = 0;

    //! Gets the mass of the grid for this evaluator
    virtual T getMass() const = 0;

    //! Gets the cell softening scale of the grid for this evaluator
    virtual T getEps() const = 0;

    //! Returns the particle at index id on the grid, without offset
    virtual Particle <T> getParticleNoOffset(size_t id) const = 0;

    //! Returns a reference to the grid for this evaluator
    const grids::Grid<T> &getGrid() const {
      return grid;
    }

    //! Returns the particle at index id on the grid for this evaluator
    Particle <T> getParticle(size_t id) const {
      Particle<T> particle = getParticleNoWrap(id);
      particle.pos = grid.wrapPoint(particle.pos);
      return particle;
    }
  };


  /*! \class ParticleGenerator
      \brief A class to generate particles starting from cells of grids.
       Note this is a pure virtual class because we require a specific implementation
       (currently Zeldovich is the only option!)
   */
  template<typename GT>
  class ParticleGenerator : public std::enable_shared_from_this<ParticleGenerator<GT>> {
  public:
    using T = tools::datatypes::strip_complex<GT>;

  protected:
    const grids::Grid<T> &grid; //!< Grid to generate particles for

    //! Constructor, requires only the grid
    ParticleGenerator(const grids::Grid<T> &grid) : grid(grid) {

    }

  public:
    //! Returns the grid for this generator
    const grids::Grid<T> &getGrid() const {
      return grid;
    }

    //! Recalculates the position and velocity offsets for particles on this grid
    virtual void recalculate() = 0;

    //! Returns a vector of the fields required to generate position and velocity offsets
    virtual std::vector<std::shared_ptr<fields::Field<GT>>> getGeneratedFields() = 0;

  };


}
#endif
