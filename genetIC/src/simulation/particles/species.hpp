#ifndef _SPECIES_HPP
#define _SPECIES_HPP

#include <map>
#include <memory>
//#include "multilevelgenerator.hpp"

namespace particle {
  template<typename T>
  class AbstractMultiLevelParticleGenerator;

  enum species {
    dm = 0, baryon = 1, whitenoise = 2, unknown = 3, all = 4
  };

  std::istream &operator>>(std::istream &inputStream, species &sp) {
    std::string s;
    inputStream >> s;
    try {
      int i = std::stoi(s);
      sp = static_cast<species>(i);
    } catch (std::invalid_argument e) {
      if (s == "dm") {
        sp = species::dm;
      } else if (s == "baryon" || s == "gas") {
        sp = species::baryon;
      } else {
        inputStream.setstate(std::ios::failbit);
      }
    }
    return inputStream;
  }


  std::ostream &operator<<(std::ostream &outputStream, const species &sp) {
    switch (sp) {
      case species::dm:
        outputStream << "dm";
        break;
      case species::baryon:
        outputStream << "baryon";
        break;
      case species::whitenoise:
        outputStream << "whitenoise";
        break;
      default:
        outputStream << "unknown";
    }
    return outputStream;
  }

  template<typename GridDataType>
  using SpeciesToGeneratorMap = std::map<species, std::shared_ptr<AbstractMultiLevelParticleGenerator<GridDataType>>>;

}

#endif