#ifndef IC_MODIFICATION_HPP
#define IC_MODIFICATION_HPP

#include <src/tools/data_types/complex.hpp>
#include <memory>

#include <src/simulation/field/multilevelfield.hpp>
#include <src/cosmology/parameters.hpp>

namespace modifications{

	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class Modification : public std::enable_shared_from_this<Modification<DataType,T>>{
	private:
		T target;

	protected:
		multilevelcontext::MultiLevelContextInformation<DataType> &underlying;
		const cosmology::CosmologicalParameters<T> &cosmology;


	public:
		Modification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
								 const cosmology::CosmologicalParameters<T> &cosmology_): underlying(underlying_), cosmology(cosmology_){};

		virtual T calculateCurrentValue(fields::MultiLevelField<DataType>* /* field */,
																		multilevelcontext::MultiLevelContextInformation<DataType>& /*&underlying*/) = 0;

		T getTarget(){
			return target;
		}

		void setTarget(T target_){
			target = target_;
		}

	};

};




#endif
