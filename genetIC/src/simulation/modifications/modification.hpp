#ifndef IC_MODIFICATION_HPP
#define IC_MODIFICATION_HPP

#include <src/tools/data_types/complex.hpp>


namespace modifications{

	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class Modification {
	private:
		T target;

	public:
		virtual T calculateCurrentValue(fields::MultiLevelField<DataType> /* field */) = 0;

		T getTarget(){
			return target;
		}

		void setTarget(T target_){
			target = target_;
		}

	};

};




#endif
