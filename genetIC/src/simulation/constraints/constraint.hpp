#ifndef IC_CONSTRAINT_HPP
#define IC_CONSTRAINT_HPP

#include <string>
#include <src/tools/data_types/complex.hpp>
#include <src/simulation/field/multilevelfield.hpp>

namespace constraints {

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class Constraint {

	//TODO Useless attribute ?
	private:
		int order;

	public:
		std::string name;
		std::string type;
		T target;
		T existing;

		Constraint() {
			name = "";
			type = "";
			target = 0;
			order = 0;
		}

		Constraint(std::string name_in, std::string type_in, T target_in, T exisiting_in){
			name = name_in;
			type = type_in;
			target = target_in;
			existing = exisiting_in;

			if (isQuadratic(name)){order = 2;}
			else if (isLinear(name)){order = 1;}
			else {throw std::runtime_error("Only linear and quadratic constraints are supported");}

		}

		bool isQuadratic(std::string name) {
			if (name == "variance") {
				return true;
			}
		}

		bool isLinear(std::string name) {
			if (name == "overdensity" || name == "phi" || name == "lx" || name == "ly" || name == "lz") {
				return true;
			}
		}

		int getOrder(){
			return order;
		}

	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearConstraint : public Constraint<DataType, T> {
	public:
		fields::ConstraintField<DataType>* alpha;

		LinearConstraint(std::string name_in, std::string type_in, T target_in, T existing_in):
				Constraint<DataType,T>(name_in,type_in,target_in,existing_in), alpha(NULL){
		}

	};


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class QuadraticConstraint : public Constraint<DataType, T> {
	public:
		std::vector<fields::ConstraintField<DataType>*> alphas;
		int initialNumberSteps;
		float filter_scale;
		float precision;

		QuadraticConstraint(std::string name_in, std::string type_in, T target_in,
												T exisitng_in, float precision_in, float filter_in, int number_in):
				Constraint<DataType,T>(name_in,type_in,target_in,existing_in){
				filter_scale = filter_in;
				initialNumberSteps = number_in;
				precision = precision_in;
				alphas = std::vector<fields::ConstraintField<DataType>*>();
		};



	};


}


#endif
