//
// Created by Andrew Pontzen on 16/07/2016.
//

#ifndef IC_PARTICLE_HPP
#define IC_PARTICLE_HPP

#include "coordinate.hpp"

template<typename T>
class Particle {
public:
  Coordinate<T> pos;
  Coordinate<T> vel;
  T mass;
  T soft;

  Particle(int) : pos(T(0)), vel(T(0)), mass(0), soft(0) {

  }

  Particle() {}

};

#endif //IC_PARTICLE_HPP
