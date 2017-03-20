//
// Created by Andrew Pontzen on 23/12/2016.
//

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
    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForGrid(const grids::Grid<T> &grid) =0;

    const particle::ParticleGenerator<GridDataType> &getGeneratorForGrid(const grids::Grid<T> &grid) const {
      return const_cast<AbstractMultiLevelParticleGenerator<GridDataType>*>(this)->getGeneratorForGrid(grid);
    }
  };

  template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
  class NullMultiLevelParticleGenerator: public AbstractMultiLevelParticleGenerator<GridDataType> {
  public:

    NullMultiLevelParticleGenerator() {

    }

    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t level) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }

    virtual particle::ParticleGenerator<GridDataType> &getGeneratorForGrid(const grids::Grid<T> &grid) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }
  };

  template<typename A,typename B,typename C>
  class MultiLevelParticleGenerator;



  template<typename GridDataType, typename T>
  void initialiseParticleGeneratorBasedOnTemplate(MultiLevelParticleGenerator<GridDataType, ZeldovichParticleGenerator<GridDataType>, T> &generator) {
    using ZPG=ZeldovichParticleGenerator<GridDataType>;
    size_t nlevels = generator.context.getNumLevels();

    generator.outputField.toFourier();

    if (nlevels == 0) {
      throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
    } else if (nlevels == 1) {
      generator.pGenerators.emplace_back(std::make_shared<ZPG>(generator.outputField.getFieldForLevel(0), generator.cosmoParams));
    } else if (nlevels >= 2) {
      cerr << "Zeldovich approximation on successive levels...";

#define NEWMETHOD

#ifndef NEWMETHOD
      // remove all the unwanted k-modes from grids, but keep a record of them so we can reinstate them on the
      // coarser levels later
      auto residuals = outputField.getHighKResiduals();
      outputField.applyFilters();
#endif

      for(size_t level=0; level<nlevels; ++level)
        generator.pGenerators.emplace_back(std::make_shared<ZPG>(generator.outputField.getFieldForLevel(level), generator.cosmoParams));

      cerr << "Interpolating low-frequency information into zoom regions...";

      for(size_t level=1; level<nlevels; ++level) {

#ifdef NEWMETHOD
        // remove the low-frequency information from this level
        generator.outputField.getFieldForLevel(level).applyFilter(generator.outputField.getHighPassFilterForLevel(level));

        // replace with the low-frequency information from the level below

        generator.outputField.getFieldForLevel(level).addFieldFromDifferentGridWithFilter(generator.outputField.getFieldForLevel(level - 1),
                                                                                          generator.outputField.getLowPassFilterForLevel(level-1));


        generator.pGenerators[level]->applyFilter(generator.outputField.getHighPassFilterForLevel(level));


        generator.pGenerators[level]->addFieldFromDifferentGridWithFilter(*generator.pGenerators[level - 1],
                                                                          generator.outputField.getLowPassFilterForLevel(level-1));


#else
        generator.outputField.getFieldForLevel(level).addFieldFromDifferentGrid(generator.outputField.getFieldForLevel(level - 1));
        generator.pGenerators[level]->addFieldFromDifferentGrid(*generator.pGenerators[level - 1]);
#endif


      }

      generator.outputField.toFourier();

#ifdef NEWMETHOD
      generator.outputField.setStateRecombined();
#else
      cerr << "Re-introducing high-k modes into low-res region...";

      generator.outputField.recombineHighKResiduals(residuals);

      for(size_t level=0; level<nlevels-1; ++level) {
        generator.pGenerators[level]->recalculate();
      }

#endif

      for(size_t i=0; i<nlevels; ++i)
        generator.pGenerators[i]->toReal();

      cerr << "done." << endl;

    }
  }



  template<typename GridDataType, typename TParticleGenerator, typename T=tools::datatypes::strip_complex<GridDataType>>
  class MultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<GridDataType>
  {
  protected:
    fields::OutputField<GridDataType> & outputField;
    const multilevelcontext::MultiLevelContextInformation<GridDataType> &context;
    std::vector<std::shared_ptr<TParticleGenerator>> pGenerators;
    const cosmology::CosmologicalParameters<T> &cosmoParams;


    void initialise() {
      // Dispatch to correct initialiser.
      // This is required because partial specialization is not allowed at class level.
      initialiseParticleGeneratorBasedOnTemplate(*this);
    }


    friend void initialiseParticleGeneratorBasedOnTemplate<>(MultiLevelParticleGenerator<GridDataType, TParticleGenerator, T> &generator);




  public:

    MultiLevelParticleGenerator(fields::OutputField<GridDataType> &field, const cosmology::CosmologicalParameters<T> &params) :
      outputField(field),
      context(field.getContext()),
      cosmoParams(params) {
      initialise();
    }

    particle::ParticleGenerator<GridDataType> &getGeneratorForLevel(size_t level) override  {
      return *(pGenerators[level]);
    }

    particle::ParticleGenerator<GridDataType> &getGeneratorForGrid(const grids::Grid<T> &grid) override {
      for(size_t i=0; i<outputField.getNumLevels(); i++) {
        if(grid.pointsToGrid(&(outputField.getFieldForLevel(i).getGrid())))
          return getGeneratorForLevel(i);
      }
      throw std::runtime_error("Attempt to get a generator for a grid which is not related to the field");
    }


  };






}



#endif //IC_MULTILEVELGENERATOR_HPP_HPP
