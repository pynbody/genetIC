#ifndef _SPECIES_HPP
#define _SPECIES_HPP

#include <map>
#include <memory>
//#include "multilevelgenerator.hpp"

namespace particle {
  template<typename T>
  class AbstractMultiLevelParticleGenerator;

  enum species {
    dm=0, baryon=1
  };

  template<typename GridDataType>
  using SpeciesToGeneratorMap = std::map<species, std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>>>;

}

#endif