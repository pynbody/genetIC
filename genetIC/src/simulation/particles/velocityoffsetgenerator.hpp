//
// Created by Andrew Pontzen on 01/10/2018.
//

#ifndef IC_VELOCITYOFFSETGENERATOR_HPP
#define IC_VELOCITYOFFSETGENERATOR_HPP

#include "multilevelgenerator.hpp"
#include "generator.hpp"
#include <memory>


namespace particle {
    /*! \class VelocityOffsetEvaluator
        \brief Class to add velocity offsets to generated particles. Requires an underlying evaluator to actually work
    */
    template<typename GridDataType>
    class VelocityOffsetEvaluator : public ParticleEvaluator<GridDataType> {
    protected:
        std::shared_ptr<ParticleEvaluator<GridDataType>> underlying; //!< Evaluator to be used in generating the particles in question
        Coordinate<GridDataType> velOffset; //!< Velocity offset to be added to all particles

    public:
        /*! \brief Constructor from the underlying evaluator and required velocity offset
            \param underlying - underlying particle evaluator
            \param velOffset - velocity offset to be added to all particles
        */
        VelocityOffsetEvaluator(std::shared_ptr<ParticleEvaluator<GridDataType>> underlying,
                Coordinate<GridDataType> velOffset) : ParticleEvaluator<GridDataType>(underlying->getGrid()),
                        underlying(underlying), velOffset(velOffset) { }

        //! Gets the particle evaluated by the underlying evaluator, and adds a velocity offset to it, without wrapping
        Particle<GridDataType> getParticleNoWrap(size_t id) const override {
          auto output = underlying->getParticleNoWrap(id);
          output.vel+=velOffset;
          return output;
        }

        //! Gets the particle evaluated by the underlying evaluator, and adds a velocity offset to it, without offsetting position
        Particle<GridDataType> getParticleNoOffset(size_t id) const override {
          auto output = underlying->getParticleNoOffset(id);
          output.vel+=velOffset;
          return output;
        }

        //! Gets the mass of the particle
        GridDataType getMass() const override {
          return underlying->getMass();
        }

        //! Gets the cell softening scale of the particle
        GridDataType getEps() const override {
          return underlying->getEps();
        }
    };

    /*! \class VelocityOffsetMultiLevelParticleGenerator
        \brief Generator that creates particles with a velocity offset added to the particle velocities

        Note that this requires an underlying generator to work - it just adds on the velocity at the end.
    */
    template<typename GridDataType>
    class VelocityOffsetMultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
    protected:
        std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>> underlying; //!< Underlying particle generator
        Coordinate<GridDataType> velOffset; //!< Velocity offset to be added to all particles
        using EvaluatorType = VelocityOffsetEvaluator<GridDataType>;

    public:
        /*! \brief Constructor, using the underlying generator, and specified velocity offset
            \param underlying - required particle generator
            \param velOffset - vector to offset velocities by
        */
        VelocityOffsetMultiLevelParticleGenerator( std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>> underlying,
                                                   Coordinate<GridDataType> velOffset) : underlying(underlying), velOffset(velOffset)
        {

        }

        //! Returns the underlying particle generator for the specified level
        ParticleGenerator <GridDataType> &getGeneratorForLevel(size_t level) override {
          return underlying->getGeneratorForLevel(level);
        }

        //! Returns a velocity evaluator which actually performs the velocity offsetting.
        std::shared_ptr<ParticleEvaluator<GridDataType>>
        makeEvaluatorForGrid(const grids::Grid<GridDataType> &grid) override {
          return std::make_shared<EvaluatorType>(underlying->makeEvaluatorForGrid(grid), velOffset);
        }

    };
}

#endif //IC_VELOCITYOFFSETGENERATOR_HPP
