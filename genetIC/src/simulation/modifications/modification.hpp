#ifndef IC_MODIFICATION_HPP
#define IC_MODIFICATION_HPP

#include <src/tools/data_types/complex.hpp>
#include <memory>


namespace modifications{

	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class Modification : public std::enable_shared_from_this<Modification<DataType,T>>{
	private:
		T target;

	public:
		virtual T calculateCurrentValue(fields::MultiLevelField<DataType>* /* field */,
																		multilevelcontext::MultiLevelContextInformation<DataType>& /*&underlying*/,
																		cosmology::CosmologicalParameters<T>& /* &cosmology*/) = 0;

		T getTarget(){
			return target;
		}

		void setTarget(T target_){
			target = target_;
		}

	};

};




#endif
