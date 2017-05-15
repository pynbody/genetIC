#ifndef IC_CONSTRAINTMANAGER_HPP
#define IC_CONSTRAINTMANAGER_HPP

#include <string>
#include "multilevelconstraintgenerator.hpp"
#include "constraintapplicator.hpp"
#include "constraint.hpp"

namespace constraints {

	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class ConstraintManager {

	public:
		MultiLevelConstraintGenerator<DataType,T> generator;
		ConstraintApplicator<DataType,T> applicator;
		fields::OutputField<DataType> *outputField;
		std::vector<LinearConstraint<DataType,T>> linearList;
		std::vector<QuadraticConstraint<DataType,T>> quadraticList;



		ConstraintManager(multilevelcontext::MultiLevelContextInformation<DataType> multiLevelContext,
											cosmology::CosmologicalParameters<T> cosmology,
											fields::OutputField<DataType> *outputField_in):
											applicator(&multiLevelContext, outputField_in),
											generator(multiLevelContext, cosmology){
											outputField = outputField_in;
		}




		void addConstrainToLinearList(std::string name, std::string type, float target){

			bool relative = isRelative(type);

			if (name != "overdensity" || name != "phi" || name != "lx" || name != "ly" || name != "lz"){
				throw std::runtime_error(name + "is not an implemented linear constraint'");
			}

			auto existing = calculateCurrentValue(name, *outputField);
			if (relative) target *= existing;
			LinearConstraint<DataType,T> constraint = LinearConstraint<DataType,T>(name, type, target, existing);
			auto covector = generator.calcConstraintForAllLevels(name);
			constraint.alpha = &covector;
			linearList.push_back(constraint);
		}



		void addConstrainToQuadList(std::string name, std::string type, float target,
																int initNumberSteps, float precision, float filterscale){

			bool relative = isRelative(type);

			if (name != "variance"){
				throw std::runtime_error(name + "is not an implemented quadratic constraint'");
			}

			auto existing = calculateCurrentValue(name,outputField);
			if (relative) target *= existing;
			QuadraticConstraint<DataType,T> constraint = QuadraticConstraint<DataType,T>(name, type, target, existing,
																																									 initNumberSteps, precision, filterscale);
		}

		T calculateCurrentValue(std::string name, fields::OutputField<DataType> field){
			auto covector = generator.calcConstraintForAllLevels(name);
			T value = covector.innerProduct(field);
			return value;
		}

		T calculateCurrentValue(LinearConstraint<DataType,T> constraint, fields::OutputField<DataType> field){
			auto covector = generator.calcConstraintForAllLevels(constraint.name);
			constraint.alpha = &covector;

		}

		void applyAllConstraints(){
			for (auto it = linearList.begin(); it != linearList.end(); ++it){
				applicator.add_constraint(linearList[*it].alpha, linearList[*it].target, linearList[*it].existing);
			}
			applicator.applyConstraints();
		}


	private:
		bool isRelative(std::string type){
			bool relative = false;
			if (strcasecmp(type.c_str(), "relative") == 0) {
				relative = true;
			} else if (strcasecmp(type.c_str(), "absolute") != 0) {
				throw std::runtime_error("Constraint type must be either 'relative' or 'absolute'");
			}
			return relative;
		}



	};
}



#endif
