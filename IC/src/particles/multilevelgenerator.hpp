//
// Created by Andrew Pontzen on 23/12/2016.
//

#ifndef IC_MULTILEVELGENERATOR_HPP_HPP
#define IC_MULTILEVELGENERATOR_HPP_HPP

#include "../multilevelfield.hpp"
#include "../cosmo.hpp"
#include "generator.hpp"
#include "zeldovich.hpp"

namespace particle {
  template<typename T>
  class AbstractMultiLevelParticleGenerator :
    public std::enable_shared_from_this<AbstractMultiLevelParticleGenerator<T>> {

  public:
    virtual particle::ParticleGenerator<T> &getGeneratorForLevel(size_t level) = 0;
    virtual particle::ParticleGenerator<T> &getGeneratorForGrid(const Grid<T> &grid) =0;

    const particle::ParticleGenerator<T> &getGeneratorForGrid(const Grid<T> &grid) const {
      return const_cast<AbstractMultiLevelParticleGenerator<T>*>(this)->getGeneratorForGrid(grid);
    }
  };

  template<typename T>
  class NullMultiLevelParticleGenerator: public AbstractMultiLevelParticleGenerator<T> {
  public:

    NullMultiLevelParticleGenerator() {

    }

    virtual particle::ParticleGenerator<T> &getGeneratorForLevel(size_t level) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }

    virtual particle::ParticleGenerator<T> &getGeneratorForGrid(const Grid<T> &grid) override {
      throw std::runtime_error("Attempt to generate particles before they have been calculated");
    }
  };

  template<typename T, typename TParticleGenerator>
  class MultiLevelParticleGenerator : public AbstractMultiLevelParticleGenerator<T>
  {
  protected:
    OutputField<std::complex<T>> & outputField;
    const MultiLevelContextInformation<T> &context;
    std::vector<std::shared_ptr<TParticleGenerator>> pGenerators;
    const CosmologicalParameters<T> &cosmoParams;

    void initialise() {
      // TODO: add specialisations
      assert(false);
    }


  public:

    MultiLevelParticleGenerator(OutputField<std::complex<T>> &field, const CosmologicalParameters<T> &params) :
      outputField(field),
      context(field.getContext()),
      cosmoParams(params) {
      initialise();
    }

    particle::ParticleGenerator<T> &getGeneratorForLevel(size_t level) override  {
      return *(pGenerators[level]);
    }

    particle::ParticleGenerator<T> &getGeneratorForGrid(const Grid<T> &grid) override {
      for(size_t i=0; i<outputField.getNumLevels(); i++) {
        if(grid.pointsToGrid(&(outputField.getFieldForLevel(i).getGrid())))
          return getGeneratorForLevel(i);
      }
      throw std::runtime_error("Attempt to get a generator for a grid which is not related to the field");
    }


  };


  // TODO: why can't I do partial specialization here?
  template<>
  void MultiLevelParticleGenerator<double, ZeldovichParticleGenerator<double>>::initialise() {
    using ZPG=ZeldovichParticleGenerator<double>;
    size_t nlevels = context.getNumLevels();

    outputField.toFourier();

    if (nlevels == 0) {
      throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
    } else if (nlevels == 1) {
      pGenerators.emplace_back(std::make_shared<ZPG>(outputField.getFieldForLevel(0), cosmoParams));
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
        pGenerators.emplace_back(std::make_shared<ZPG>(outputField.getFieldForLevel(level), cosmoParams));

      cerr << "Interpolating low-frequency information into zoom regions...";

      for(size_t level=1; level<nlevels; ++level) {

#ifdef NEWMETHOD
        // remove the low-frequency information from this level
        outputField.getFieldForLevel(level).applyFilter(outputField.getHighPassFilterForLevel(level));

        // replace with the low-frequency information from the level below

        outputField.getFieldForLevel(level).addFieldFromDifferentGridWithFilter(outputField.getFieldForLevel(level - 1),
                                                                                outputField.getLowPassFilterForLevel(level-1));


        pGenerators[level]->applyFilter(outputField.getHighPassFilterForLevel(level));


        pGenerators[level]->addFieldFromDifferentGridWithFilter(*pGenerators[level - 1],
                                                                outputField.getLowPassFilterForLevel(level-1));


#else
        outputField.getFieldForLevel(level).addFieldFromDifferentGrid(outputField.getFieldForLevel(level - 1));
        pGenerators[level]->addFieldFromDifferentGrid(*pGenerators[level - 1]);
#endif


      }

      outputField.toFourier();

#ifdef NEWMETHOD
      outputField.setStateRecombined();
#else
      cerr << "Re-introducing high-k modes into low-res region...";

      outputField.recombineHighKResiduals(residuals);

      for(size_t level=0; level<nlevels-1; ++level) {
        pGenerators[level]->recalculate();
      }

#endif

      for(size_t i=0; i<nlevels; ++i)
        pGenerators[i]->toReal();

      cerr << "done." << endl;

    }
  }


}



#endif //IC_MULTILEVELGENERATOR_HPP_HPP
