//
// Created by Andrew Pontzen on 23/12/2016.
//

#ifndef IC_MULTILEVELGENERATOR_HPP_HPP
#define IC_MULTILEVELGENERATOR_HPP_HPP

#include "../multilevelfield.hpp"
#include "../cosmo.hpp"
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

    outputField.toFourier();

    if (context.getNumLevels() == 0) {
      throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
    } else if (context.getNumLevels() == 1) {
      pGenerators.emplace_back(std::make_shared<ZPG>(outputField.getFieldForLevel(0), cosmoParams));
    } else if (context.getNumLevels() == 2) {
      auto residuals = outputField.getHighKResiduals();
      cerr << "Zeldovich approximation on successive levels...";
      outputField.applyFilters();
      pGenerators.emplace_back(std::make_shared<ZPG>(outputField.getFieldForLevel(0), cosmoParams));
      pGenerators.emplace_back(std::make_shared<ZPG>(outputField.getFieldForLevel(1), cosmoParams));
      cerr << "done." << endl;

      cerr << "Interpolating low-frequency information into zoom region...";
      outputField.toReal();
      outputField.getFieldForLevel(1).addFieldFromDifferentGrid(outputField.getFieldForLevel(0));
      pGenerators[1]->addFieldFromDifferentGrid(*pGenerators[0]);
      outputField.toFourier();
      // TODO: also the velocity offsets!
      cerr << "done." << endl;

      cerr << "Re-introducing high-k modes into low-res region...";
      outputField.recombineHighKResiduals(residuals);
      pGenerators[0]->recalculate();
      cerr << "done." << endl;

    } else {
      throw std::runtime_error("Too many levels!"); // TODO: general multi-level zeldovich algorithm
    }
  }


}



#endif //IC_MULTILEVELGENERATOR_HPP_HPP
