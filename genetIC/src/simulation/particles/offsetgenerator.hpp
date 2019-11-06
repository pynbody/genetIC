//
// Created by Andrew Pontzen on 01/10/2018.
//

#ifndef IC_VELOCITYOFFSETGENERATOR_HPP
#define IC_VELOCITYOFFSETGENERATOR_HPP

#include "multilevelgenerator.hpp"
#include "generator.hpp"
#include <memory>


namespace particle {
  /*! \class OffsetEvaluator
      \brief Class to add velocity offsets to generated particles from an underlying evaluator.
  */
  template<typename GridDataType>
  class OffsetEvaluator : public ParticleEvaluator<GridDataType> {
  protected:
    std::shared_ptr<
      ParticleEvaluator < GridDataType>> underlying; //!< Evaluator to be used in generating the particles in question
    Coordinate<GridDataType> velOffset; //!< Velocity offset to be added to all particles
    Coordinate<GridDataType> posOffset; //!< Position offset to be added to all particles

  public:
    /*! \brief Constructor from the underlying evaluator and required velocity offset
        \param underlying - underlying particle evaluator
        \param velOffset - velocity offset to be added to all particles
    */
    OffsetEvaluator(std::shared_ptr<ParticleEvaluator<GridDataType>> underlying,
                    Coordinate<GridDataType> posOffset, Coordinate<GridDataType> velOffset) :
                      ParticleEvaluator<GridDataType>(underlying->getGrid()),
                      underlying (underlying),
                      velOffset(velOffset), posOffset(posOffset) {}

    Particle <GridDataType> getParticleNoWrap(size_t id) const override {
      auto output = underlying->getParticleNoWrap(id);
      output.vel += velOffset;
      output.pos += posOffset;
      return output;
    }

    Particle <GridDataType> getParticleNoOffset(size_t id) const override {
      auto output = underlying->getParticleNoOffset(id);
      output.vel += velOffset;
      return output;
    }

    GridDataType getMass() const override {
      return underlying->getMass();
    }

    GridDataType getEps() const override {
      return underlying->getEps();
    }
  };

  /*! \class OffsetMultiLevelParticleGenerator
      \brief Generator that creates particles with a velocity offset added to the particle velocities

      Note that this requires an underlying generator to work - it just adds on the velocity at the end.
  */
  template<typename GridDataType>
  class OffsetMultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
  protected:
    std::shared_ptr<AbstractMultiLevelParticleGenerator < GridDataType>> underlying; //!< Underlying particle generator
    Coordinate<GridDataType> velOffset; //!< Velocity offset to be added to all particles
    Coordinate<GridDataType> posOffset; //!< Position offset to be added to all particles
    using EvaluatorType = OffsetEvaluator<GridDataType>;
    using T = tools::datatypes::strip_complex<GridDataType>;

  public:
    /*! \brief Constructor, using the underlying generator, and specified velocity offset
        \param underlying - required particle generator
        \param velOffset - vector to offset velocities by
    */
    OffsetMultiLevelParticleGenerator(std::shared_ptr<AbstractMultiLevelParticleGenerator < GridDataType>>

    underlying,
    Coordinate<GridDataType> posOffset, Coordinate<GridDataType>
    velOffset)
    :

    underlying (underlying), velOffset(velOffset), posOffset(posOffset) {

    }

    //! Returns the underlying particle generator for the specified level
    ParticleGenerator <GridDataType> &getGeneratorForLevel(size_t level) override {
      return underlying->getGeneratorForLevel(level);
    }

    //! Returns a velocity evaluator which actually performs the velocity offsetting.
    std::shared_ptr<ParticleEvaluator < GridDataType>>
    makeParticleEvaluatorForGrid(
    const grids::Grid<GridDataType> &grid
    ) override {
      return std::make_shared<EvaluatorType>(underlying->makeParticleEvaluatorForGrid(grid), posOffset, velOffset);
    }

    std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>
    makeOverdensityEvaluatorForGrid(const grids::Grid<T> &grid) override {
      return underlying->makeOverdensityEvaluatorForGrid(grid);
    }

  };
}

#endif //IC_VELOCITYOFFSETGENERATOR_HPP
