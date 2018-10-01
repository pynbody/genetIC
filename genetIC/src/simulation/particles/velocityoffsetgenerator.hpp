//
// Created by Andrew Pontzen on 01/10/2018.
//

#ifndef IC_VELOCITYOFFSETGENERATOR_HPP
#define IC_VELOCITYOFFSETGENERATOR_HPP

#include "multilevelgenerator.hpp"
#include "generator.hpp"
#include <memory>


namespace particle {
    template<typename GridDataType>
    class VelocityOffsetEvaluator : public ParticleEvaluator<GridDataType> {
    protected:
        std::shared_ptr<ParticleEvaluator<GridDataType>> underlying;
        Coordinate<GridDataType> velOffset;

    public:
        VelocityOffsetEvaluator(std::shared_ptr<ParticleEvaluator<GridDataType>> underlying,
                Coordinate<GridDataType> velOffset) : ParticleEvaluator<GridDataType>(underlying->getGrid()),
                        underlying(underlying), velOffset(velOffset) { }

        Particle<GridDataType> getParticleNoWrap(size_t id) const override {
          auto output = underlying->getParticleNoWrap(id);
          output.vel+=velOffset;
          return output;
        }

        Particle<GridDataType> getParticleNoOffset(size_t id) const override {
          auto output = underlying->getParticleNoOffset(id);
          output.vel+=velOffset;
          return output;
        }

        GridDataType getMass() const override {
          return underlying->getMass();
        }

        GridDataType getEps() const override {
          return underlying->getEps();
        }
    };

    template<typename GridDataType>
    class VelocityOffsetMultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
    protected:
        std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>> underlying;
        Coordinate<GridDataType> velOffset;
        using EvaluatorType = VelocityOffsetEvaluator<GridDataType>;

    public:
        VelocityOffsetMultiLevelParticleGenerator( std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>> underlying,
                                                   Coordinate<GridDataType> velOffset) : underlying(underlying), velOffset(velOffset)
        {

        }

        ParticleGenerator <GridDataType> &getGeneratorForLevel(size_t level) override {
          return underlying->getGeneratorForLevel(level);
        }

        std::shared_ptr<ParticleEvaluator<GridDataType>>
        makeEvaluatorForGrid(const grids::Grid<GridDataType> &grid) override {
          return std::make_shared<EvaluatorType>(underlying->makeEvaluatorForGrid(grid), velOffset);
        }

    };
}

#endif //IC_VELOCITYOFFSETGENERATOR_HPP
