//
// Created by Andrew Pontzen on 05/08/15.
//

#ifndef IC_MULTILEVELCONSTRAINTGENERATOR_H
#define IC_MULTILEVELCONSTRAINTGENERATOR_H

#include <complex.h>
#include <vector>

#include "multilevelcontext.hpp"
#include "grid.hpp"
#include "sparse.hpp"
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

  vector<complex<T>> calcConstraintVector(string name_in, int level) {
    auto ar = fieldManager.createEmptyFieldForLevel(level);

    calcConstraint(name_in, fieldManager.getGridForLevel(level), cosmology, ar);

    if (level != 0)
      inPlaceMultiply(ar, pow(fieldManager.getGridForLevel(level).dx / fieldManager.getGridForLevel(0).dx, -3.0));

    return ar;
  }

  vector<vector<complex<T>>> calcConstraintForEachLevel(string name_in, bool kspace = true) {
    auto highResConstraint = calcConstraintVector(name_in, fieldManager.getNumLevels() - 1);

    auto dataOnLevels = fieldManager.generateMultilevelFromHighResField(std::move(highResConstraint));


    if (!kspace) {
      for (auto levelData : dataOnLevels) {
        fft(levelData, levelData, -1);
      }
    }

    return dataOnLevels;
  }

  auto calcConstraintForAllLevels(string name, bool kspace = true) {
    auto dataOnLevels = calcConstraintForEachLevel(name, kspace);
    return ConstraintField<std::complex<T>>(fieldManager,std::move(dataOnLevels));
  }
};


#endif //IC_MULTILEVELCONSTRAINTGENERATOR_H
