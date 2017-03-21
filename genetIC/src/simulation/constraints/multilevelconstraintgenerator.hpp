#ifndef IC_MULTILEVELCONSTRAINTGENERATOR_H
#define IC_MULTILEVELCONSTRAINTGENERATOR_H

#include <complex.h>
#include <vector>

#include "src/simulation/multilevelcontext/multilevelcontext.hpp"
#include "src/simulation/grid/grid.hpp"
#include "src/simulation/field/field.hpp"
#include "src/tools/numerics/vectormath.hpp"
#include "src/simulation/constraints/onelevelconstraintgenerator.hpp"
#include "src/simulation/field/multilevelfield.hpp"



namespace constraints {

  /** This class takes responsibility for calculating constraint covectors across all levels that a field is
    * defined on. */
  template<typename DataType, typename T = tools::datatypes::strip_complex<DataType>>
  class MultiLevelConstraintGenerator {


  protected:
    multilevelcontext::MultiLevelContextInformation<DataType> &fieldManager;
    cosmology::CosmologicalParameters<T> &cosmology;

    fields::Field<DataType> calcConstraintCovectorNormalisedForLevel(std::string name_in, int level) {
      /* Get a named constraint covector on a particular level. Normalise such that taking the multi-level inner product
       * will return the correct value (the same as it would if the entire box consisted just of the specified level) */

      using tools::numerics::operator/=;

      // calculate on the finest field
      auto returnField = calcConstraint<DataType>(name_in, fieldManager.getGridForLevel(level), cosmology);

      // When taking inner products, cells on fine levels get a weight proportional to their volume element. However,
      // we want to regard the covector just calculated on the finest level as definitive, so this weight
      // must be pre-divided out.
      if (level != 0)
        returnField.getDataVector() /= fieldManager.getWeightForLevel(level);

      return std::move(returnField);
    }

  public:
    MultiLevelConstraintGenerator(multilevelcontext::MultiLevelContextInformation<DataType> &fieldManager,
                                  cosmology::CosmologicalParameters<T> &cosmology) :
      fieldManager(fieldManager), cosmology(cosmology) {

    }


    auto calcConstraintForAllLevels(std::string name) {
      /* Calculate the multi-level covector for a named constraint */
      auto highResConstraint = calcConstraintCovectorNormalisedForLevel(name, fieldManager.getNumLevels() - 1);
      return fieldManager.generateMultilevelFromHighResField(std::move(highResConstraint));
    }
  };
}

#endif //IC_MULTILEVELCONSTRAINTGENERATOR_H
