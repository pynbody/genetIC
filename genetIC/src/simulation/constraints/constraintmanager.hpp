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
		fields::OutputField<DataType>* outputField;
		MultiLevelConstraintGenerator<DataType,T> generator;
		ConstraintApplicator<DataType,T> applicator;

		std::vector<LinearConstraint<DataType,T>> linearList;
		std::vector<QuadraticConstraint<DataType,T>> quadraticList;



		ConstraintManager(multilevelcontext::MultiLevelContextInformation<DataType>* multiLevelContext,
											cosmology::CosmologicalParameters<T> cosmology,
											fields::OutputField<DataType>* outputField_in):
											outputField(outputField_in),
											generator(*multiLevelContext, cosmology),
											applicator(multiLevelContext, outputField) {}


		T calculateCurrentValueByName(std::string name){
			auto covector = generator.calcConstraintForAllLevels(name);
			covector.toFourier();
			T value = covector.innerProduct(*outputField).real();
			return value;
		}

		void addConstrainToLinearList(std::string name, std::string type, float target){

			bool relative = isRelative(type);

			if (name == "overdensity" || name == "phi" || name == "lx" || name == "ly" || name == "lz") {

				auto existing = calculateCurrentValueByName(name);
				if (relative) target *= existing;
				LinearConstraint<DataType, T> constraint = LinearConstraint<DataType, T>(name, type, target, existing);
				auto covector = generator.calcConstraintForAllLevels(name);
				covector.toFourier();
				constraint.setAlpha(&covector);
				linearList.push_back(std::move(constraint));
			} else{
				throw std::runtime_error(name + "" +  "is not an implemented linear constraint'");
			}
		}



		void addConstrainToQuadList(std::string name, std::string type, float target,
																int initNumberSteps, float precision, float filterscale){

			bool relative = isRelative(type);

			if (name != "variance"){
				throw std::runtime_error(name + "is not an implemented quadratic constraint'");
			}

			T existing = calculateCurrentValueByName(name);
			if (relative) target *= existing;
			QuadraticConstraint<DataType,T> constraint = QuadraticConstraint<DataType,T>(name, type, target, existing,
																																									 initNumberSteps, precision, filterscale);
		}

		void applyAllConstraints(){
			std::cout << "Starting loop" << std::endl;
			for (auto it = linearList.begin(); it != linearList.end(); ++it){
				auto a = (*it).alpha;
				std::cout << "Ending loop" << std::endl;
				applicator.add_constraint(std::move(*a), (*it).target, (*it).existing);
			}
			std::cout << "Entering applicator.applyConstraints" << std::endl;
			applicator.applyConstraints();
			std::cout << "Leaving applicator.applyConstraints" << std::endl;
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
