#ifndef _SPECIES_HPP
#define _SPECIES_HPP

#include <map>
#include <memory>
#include "multilevelgenerator.hpp"

namespace particle {
  enum species {
    dm, baryon
  };

  template<typename GridDataType>
  using SpeciesToGeneratorMap = std::map<species, std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>>>;

}

#endif