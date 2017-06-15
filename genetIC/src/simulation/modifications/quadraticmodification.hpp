#ifndef IC_QUADRATICMODIFICATION_HPP
#define IC_QUADRATICMODIFICATION_HPP

#include <src/simulation/modifications/modification.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class QuadraticModification : public Modification<DataType, T> {
	private:
		int initNumberSteps;
		float precision;

	public:
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class VarianceModification : public QuadraticModification<DataType, T> {
	private:
		float scale;

	public:
	};
}



#endif
