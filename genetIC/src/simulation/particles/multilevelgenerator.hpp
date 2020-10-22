#ifndef IC_MULTILEVELGENERATOR_HPP_HPP
#define IC_MULTILEVELGENERATOR_HPP_HPP

#include "src/simulation/field/multilevelfield.hpp"
#include "src/cosmology/parameters.hpp"
#include "src/simulation/particles/generator.hpp"
#include "src/simulation/particles/zeldovich.hpp"
#include "src/simulation/field/evaluator.hpp"

#include <memory>

namespace particle {
  using std::cerr;
  using std::endl;

  /*! \class AbstractMultiLevelParticleGenerator
      \brief Base class for all multi-level particle generator
  */
  template<typename GridDataType>
  class AbstractMultiLevelParticleGenerator :
    public std::enable_shared_from_this<AbstractMultiLevelParticleGenerator<GridDataType>> {

  public:
    using T = tools::datatypes::strip_complex<GridDataType>;

    //! Returns a reference to the generator for the specified level
    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t level) = 0;

    //! Creates a particle evaluator for the specified grid and returns a pointer to it
    virtual std::shared_ptr<particle::ParticleEvaluator<GridDataType>>
    makeParticleEvaluatorForGrid(const grids::Grid<T> &grid) = 0;

    virtual std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>
    makeOverdensityEvaluatorForGrid(const grids::Grid<T> &grid) = 0;

    //! Creates a particle evaluator for the specified grid and returns a constant pointer to it
    std::shared_ptr<const particle::ParticleEvaluator<GridDataType>>
    makeParticleEvaluatorForGrid(const grids::Grid<T> &grid) const {
      return const_cast<AbstractMultiLevelParticleGenerator<GridDataType> *>(this)->makeParticleEvaluatorForGrid(grid);
    }
  };

  /*! \class NullMultiLevelParticleGenerator
      \brief Class for undefined particle generators - intended to be a place-holder
  */
  template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
  class NullMultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
  public:

    //! Default constructor.
    NullMultiLevelParticleGenerator() {

    }

    //! Throws an error, because we can't generate particles before we have defined the generator
    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t /*level*/) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }

    //! Throws an error, because we can't generate particles before we have defined the generator
    virtual std::shared_ptr<particle::ParticleEvaluator<GridDataType>>
    makeParticleEvaluatorForGrid(const grids::Grid<T> &) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }

    std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>
    makeOverdensityEvaluatorForGrid(const grids::Grid<T> &grid) override {
      throw std::runtime_error("Attempt to get overdensity before it has been calculated");
    }

  };

  template<typename A, typename B, typename C>
  class MultiLevelParticleGenerator;


  //! Function to initialise Zeldovich particle generators. Separated from main class due to lack of c++ partial template specialisation
  template<typename GridDataType, typename T>
  void initialiseParticleGeneratorBasedOnTemplate(
    MultiLevelParticleGenerator<GridDataType, ZeldovichParticleGenerator<GridDataType>, T> &generator) {
    using ZPG=ZeldovichParticleGenerator<GridDataType>;
    size_t nlevels = generator.context.getNumLevels();

    logging::entry() << "Calculating particles from overdensity fields..." << endl;

    generator.overdensityField.toFourier();

    if (nlevels == 0) {
      throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
    } else if (nlevels == 1) {
      generator.pGenerators.emplace_back(
        std::make_shared<ZPG>(generator.overdensityField.getFieldForLevel(0)));
    } else if (nlevels >= 2) {


      for (size_t level = 0; level < nlevels; ++level)
        generator.pGenerators.emplace_back(
          std::make_shared<ZPG>(generator.overdensityField.getFieldForLevel(level)));

      logging::entry() << "Combining information from different levels..." << endl;

      auto filters = generator.overdensityField.getFilters();

      for (size_t level = 1; level < nlevels; ++level) {

        // remove the low-frequency information from this level
        generator.overdensityField.getFieldForLevel(level).applyFilter(
          filters.getHighPassFilterForLevel(level));
        generator.pGenerators[level]->applyFilter(filters.getHighPassFilterForLevel(level));

        // replace with the low-frequency information from the level below
        generator.overdensityField.getFieldForLevel(level).addFieldFromDifferentGridWithFilter(
          generator.overdensityField.getFieldForLevel(level - 1),
          filters.getLowPassFilterForLevel(level - 1));
        generator.pGenerators[level]->addFieldFromDifferentGridWithFilter(*generator.pGenerators[level - 1],
                                                                          filters.getLowPassFilterForLevel(level - 1));


      }

      generator.overdensityField.getContext().setLevelsAreCombined();

      for (size_t i = 0; i < nlevels; ++i)
        generator.pGenerators[i]->toReal();


    }
  }

  //! Function to to make a Zeldovich particle evaluator. Separated from main class due to lack of c++ partial template specialisation
  template<typename GridDataType, typename T>
  std::shared_ptr<particle::ParticleEvaluator<GridDataType>> makeParticleEvaluatorBasedOnTemplate(
    MultiLevelParticleGenerator<GridDataType, ZeldovichParticleGenerator<GridDataType>, T> &generator,
    const grids::Grid<T> &grid, T epsNorm = 0.01075) {
    auto fieldEvaluators = generator.getOutputFieldEvaluatorsForGrid(grid);
    return std::make_shared<ZeldovichParticleEvaluator<GridDataType>>(fieldEvaluators[0], fieldEvaluators[1],
                                                                      fieldEvaluators[2],
                                                                      grid, generator.cosmoParams, epsNorm);
  };


  /*! \class MultiLevelParticleGenerator
      \brief Main class object used to handle generator for a multi-level system
  */
  template<typename GridDataType, typename TParticleGenerator, typename T=tools::datatypes::strip_complex<GridDataType>>
  class MultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType> {
  protected:
    using ContextType = multilevelgrid::MultiLevelGrid<GridDataType>;
    fields::OutputField<GridDataType> &overdensityField; //!< Reference to the multi-level field (ie, overdensity) we need to generate particles
    const ContextType &context; //!< Reference to the multi-level context
    std::vector<std::shared_ptr<TParticleGenerator>> pGenerators; //!< Vector of generators for each level
    const cosmology::CosmologicalParameters<T> &cosmoParams; //!< Cosmological parameters
    std::vector<std::shared_ptr<fields::MultiLevelField<GridDataType>>> outputFields; //!< Fields used to define position and velocity offsets (not overdensities!)
    T epsNorm; //!< Cell softening scale pre-factor

    //! Initialise the relevant particle generator, and gather the output fields used to define position/velocity offsets
    void initialise() {
      // Dispatch to correct initialiser.
      // This is required because partial specialization is not allowed at class level.
      initialiseParticleGeneratorBasedOnTemplate(*this);
      gatherOutputFields();
    }

    //! Return the number of output fields used to define position/velocity offsets
    size_t getNumFields() {
      std::vector<std::shared_ptr<fields::Field<GridDataType>>> fieldsOnLevel0 = pGenerators[0]->getGeneratedFields();
      return fieldsOnLevel0.size();
    }

    //! Gathers the output fields used to define position and velocity offsets, storing them in the outputFields vector
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

    //! Returns evaluator for the output fields.
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
      const grids::Grid<T> &grid, T epsNorm
    );


  public:

    /*! \brief Constructor that uses the specified overdensity field and cosmological parameters
        \param field - overdensity field needed to define particles.
        \param params - cosmological parameters.
        \param epsNorm_ - pre-factor for cell softening scale
    */
    MultiLevelParticleGenerator(fields::OutputField<GridDataType> &field,
                                const cosmology::CosmologicalParameters<T> &params, T epsNorm_ = 0.01075) :
      overdensityField(field),
      context(field.getContext()),
      cosmoParams(params) {
      initialise();
      epsNorm = epsNorm_;
    }

    //! Returns a reference to the generator on the specified level
    particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t level) override {
      return *(pGenerators[level]);
    }

    //! Returns the relevant particle evaluator
    std::shared_ptr<particle::ParticleEvaluator<GridDataType>>
    makeParticleEvaluatorForGrid(const grids::Grid<T> &grid) override {
      return makeParticleEvaluatorBasedOnTemplate(*this, grid, epsNorm);
    }

    std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>
    makeOverdensityEvaluatorForGrid(const grids::Grid<T> &grid) override {
      overdensityField.toReal();
      return fields::makeEvaluator(overdensityField, grid);
    }


  };

}

#endif
