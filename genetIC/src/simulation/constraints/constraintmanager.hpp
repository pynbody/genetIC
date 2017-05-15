#ifndef IC_CONSTRAINTMANAGER_HPP
#define IC_CONSTRAINTMANAGER_HPP

#include <string>
#include "multilevelconstraintgenerator.hpp"
#include "constraintapplicator.hpp"

namespace constraints {

	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class ConstraintManager {

	public:
		MultiLevelConstraintGenerator<DataType> generator;
		ConstraintApplicator<DataType> applicator;
		fields::OutputField<DataType> *outputField;
		std::vector<LinearConstraint<T>> linearList;
		std::vector<QuadraticConstraint<T>> quadraticList;



		ConstraintManager(multilevelcontext::MultiLevelContextInformation<DataType> multiLevelContext,
											cosmology::CosmologicalParameters<T> cosmology,
											fields::OutputField<DataType> *outputField_in):
											constraintApplicator(&multiLevelContext, &outputField_in),
											constraintGenerator(multiLevelContext, cosmology){
											outputField = outputField_in;
		}




		void addConstrainToLinearList(std::string name, std::string type, float target){

			bool relative = false;
			if (strcasecmp(type.c_str(), "relative") == 0) {
				relative = true;
			} else if (strcasecmp(type.c_str(), "absolute") != 0) {
				throw std::runtime_error("Constraint type must be either 'relative' or 'absolute'");
			}

			if (name != "overdensity" && name != "phi" && name != "lx" && name != "ly" && name != "lz"){
				throw std::runtime_error(name + "is not an implemented linear constraint'");
			}

			auto existing = calculateCurrentValue(name,outputField);
			if (relative) target *= existing;
			LinearConstraint constraint = LinearConstraint(name,type,target,existing);
			constraint.setAlphas(generator.calcConstraintForAllLevels(name));

		}

		void addConstrainToQuadList(std::string name, std::string type, float target,
																int initNumberSteps, float precision, float filterscale);

		T calculateCurrentValue(std::string name, Field);

		T calculateConstraintValue(Constraint, Field);

		void applyAllConstraints();





	};
}



#endif
