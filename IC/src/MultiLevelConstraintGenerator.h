//
// Created by Andrew Pontzen on 05/08/15.
//

#ifndef IC_MULTILEVELCONSTRAINTGENERATOR_H
#define IC_MULTILEVELCONSTRAINTGENERATOR_H

#include <complex.h>
#include <vector>

#include "multilevelcontext.hpp"
#include "grid.hpp"
#include "vectormath.hpp"
#include "onelevelconstraint.hpp"
#include "multilevelfield.hpp"

using namespace std;

template<typename T>
class MultiLevelConstraintGenerator {
protected:
  MultiLevelContextInformation<T> &fieldManager;
  CosmologicalParameters<T> &cosmology;
public:
  MultiLevelConstraintGenerator(MultiLevelContextInformation<T> &fieldManager, CosmologicalParameters<T> &cosmology) :
    fieldManager(fieldManager), cosmology(cosmology) {

  }

  Field<complex<T>, T> calcConstraintVector(string name_in, int level) {
    auto ar = fieldManager.createEmptyFieldForLevel(level);

    calcConstraint(name_in, fieldManager.getGridForLevel(level), cosmology, ar);

    if (level != 0)
      ar*=pow(fieldManager.getGridForLevel(level).dx / fieldManager.getGridForLevel(0).dx, -3.0);

    return Field<complex<T>, T>(fieldManager.getGridForLevel(level), ar, true);
  }

  auto calcConstraintForAllLevels(string name) {
    auto highResConstraint = calcConstraintVector(name, fieldManager.getNumLevels() - 1);
    return fieldManager.generateMultilevelFromHighResField(std::move(highResConstraint));
  }
};


#endif //IC_MULTILEVELCONSTRAINTGENERATOR_H
