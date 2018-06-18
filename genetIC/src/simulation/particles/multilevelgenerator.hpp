#ifndef IC_MULTILEVELGENERATOR_HPP_HPP
#define IC_MULTILEVELGENERATOR_HPP_HPP

#include "src/simulation/field/multilevelfield.hpp"
#include "src/cosmology/parameters.hpp"
#include "src/simulation/particles/generator.hpp"
#include "src/simulation/particles/zeldovich.hpp"

namespace particle {
  using std::cerr;
  using std::endl;

  template<typename GridDataType>
  class AbstractMultiLevelParticleGenerator :
      public std::enable_shared_from_this<AbstractMultiLevelParticleGenerator<GridDataType>> {

  public:
    using T = tools::datatypes::strip_complex<GridDataType>;

    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t level) = 0;

    virtual std::shared_ptr<particle::ParticleEvaluator<GridDataType>>
    makeEvaluatorForGrid(const grids::Grid<T> &grid) =0;

    std::shared_ptr<const particle::ParticleEvaluator<GridDataType>>
    makeEvaluatorForGrid(const grids::Grid<T> &grid) const {
      return const_cast<AbstractMultiLevelParticleGenerator<GridDataType> *>(this)->makeEvaluatorForGrid(grid);
    }
  };

  template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
  class NullMultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
  public:

    NullMultiLevelParticleGenerator() {

    }

    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t /*level*/) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }

    virtual std::shared_ptr<particle::ParticleEvaluator<GridDataType>>
    makeEvaluatorForGrid(const grids::Grid<T> &) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }
  };

  template<typename A, typename B, typename C>
  class MultiLevelParticleGenerator;


  template<typename GridDataType, typename T>
  void initialiseParticleGeneratorBasedOnTemplate(
      MultiLevelParticleGenerator<GridDataType, ZeldovichParticleGenerator<GridDataType>, T> &generator) {
    using ZPG=ZeldovichParticleGenerator<GridDataType>;
    size_t nlevels = generator.context.getNumLevels();

    generator.outputField.toFourier();

    if (nlevels == 0) {
      throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
    } else if (nlevels == 1) {
      generator.pGenerators.emplace_back(
          std::make_shared<ZPG>(generator.outputField.getFieldForLevel(0)));
    } else if (nlevels >= 2) {
      cerr << "Zeldovich approximation on successive levels..." << endl;


      for (size_t level = 0; level < nlevels; ++level)
        generator.pGenerators.emplace_back(
            std::make_shared<ZPG>(generator.outputField.getFieldForLevel(level)));

      cerr << "Interpolating low-frequency information into zoom regions..." << endl;

      for (size_t level = 1; level < nlevels; ++level) {

        // remove the low-frequency information from this level
        generator.outputField.getFieldForLevel(level).applyFilter(
            generator.outputField.getHighPassFilterForLevel(level));

        // replace with the low-frequency information from the level below

        generator.outputField.getFieldForLevel(level).addFieldFromDifferentGridWithFilter(
            generator.outputField.getFieldForLevel(level - 1),
            generator.outputField.getLowPassFilterForLevel(level - 1));


        generator.pGenerators[level]->applyFilter(generator.outputField.getHighPassFilterForLevel(level));


        generator.pGenerators[level]->addFieldFromDifferentGridWithFilter(*generator.pGenerators[level - 1],
                                                                          generator.outputField.getLowPassFilterForLevel(
                                                                              level - 1));


      }

      generator.outputField.toFourier();

      generator.outputField.setStateRecombined();

      for (size_t i = 0; i < nlevels; ++i)
        generator.pGenerators[i]->toReal();

      cerr << "done." << endl;

    }
  }

  template<typename GridDataType, typename T>
  std::shared_ptr<particle::ParticleEvaluator<GridDataType>> makeParticleEvaluatorBasedOnTemplate(
      MultiLevelParticleGenerator<GridDataType, ZeldovichParticleGenerator<GridDataType>, T> &generator,
      const grids::Grid<T> &grid) {
    auto fieldEvaluators = generator.getOutputFieldEvaluatorsForGrid(grid);
    return std::make_shared<ZeldovichParticleEvaluator<GridDataType>>(fieldEvaluators[0], fieldEvaluators[1],
                                                                      fieldEvaluators[2],
                                                                      grid, generator.cosmoParams);
  };


  template<typename GridDataType, typename TParticleGenerator, typename T=tools::datatypes::strip_complex<GridDataType>>
  class MultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
  protected:
    using ContextType = multilevelcontext::MultiLevelContextInformation<GridDataType>;
    fields::OutputField<GridDataType> &outputField;
    const ContextType &context;
    std::vector<std::shared_ptr<TParticleGenerator>> pGenerators;
    const cosmology::CosmologicalParameters<T> &cosmoParams;
    std::vector<std::shared_ptr<fields::MultiLevelField<GridDataType>>> outputFields;


    void initialise() {
      // Dispatch to correct initialiser.
      // This is required because partial specialization is not allowed at class level.
      initialiseParticleGeneratorBasedOnTemplate(*this);
      gatherOutputFields();
    }

    size_t getNumFields() {
      std::vector<std::shared_ptr<fields::Field<GridDataType>>> fieldsOnLevel0 = pGenerators[0]->getGeneratedFields();
      return fieldsOnLevel0.size();
    }

    void gatherOutputFields() {
      size_t numFields = getNumFields();
      for (size_t field = 0; field < numFields; ++field) {
        std::vector<std::shared_ptr<fields::Field<GridDataType>>> fieldsAcrossLevels;
        for (size_t level = 0; level < context.getNumLevels(); level++) {
          fieldsAcrossLevels.push_back(pGenerators[level]->getGeneratedFields()[field]);
        }
        outputFields.emplace_back(std::make_shared<fields::MultiLevelField<GridDataType>>(
            const_cast<ContextType &>(context), fieldsAcrossLevels));
      }
    }

    std::vector<std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>>
    getOutputFieldEvaluatorsForGrid(const grids::Grid<T> &grid) {
      std::vector<std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>> fieldEvaluators;

      for (auto field : outputFields) {
        fieldEvaluators.emplace_back(fields::makeEvaluator(*field, grid));
      }

      return fieldEvaluators;
    }


    friend void initialiseParticleGeneratorBasedOnTemplate<>(
        MultiLevelParticleGenerator<GridDataType, TParticleGenerator, T> &generator);

    friend std::shared_ptr<particle::ParticleEvaluator<GridDataType>> makeParticleEvaluatorBasedOnTemplate<>(
        MultiLevelParticleGenerator<GridDataType, TParticleGenerator, T> &generator,
        const grids::Grid<T> &grid
    );


  public:

    MultiLevelParticleGenerator(fields::OutputField<GridDataType> &field,
                                const cosmology::CosmologicalParameters<T> &params) :
        outputField(field),
        context(field.getContext()),
        cosmoParams(params) {
      initialise();
    }

    particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t level) override {
      return *(pGenerators[level]);
    }

    std::shared_ptr<particle::ParticleEvaluator<GridDataType>>
    makeEvaluatorForGrid(const grids::Grid<T> &grid) override {


      return makeParticleEvaluatorBasedOnTemplate(*this, grid);
    }


  };

}

#endif
