#ifndef IC_MODIFICATIONMANAGER_HPP
#define IC_MODIFICATIONMANAGER_HPP

#include <string>
#include <src/tools/data_types/complex.hpp>
#include <src/simulation/field/multilevelfield.hpp>
#include <src/cosmology/parameters.hpp>

#include <src/simulation/modifications/linearmodification.hpp>

namespace modifications{

	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class ModificationManager {

	public:
		fields::OutputField<DataType>* outputField;
		multilevelcontext::MultiLevelContextInformation<DataType> &underlying;
		cosmology::CosmologicalParameters<T> cosmology;
		std::vector<Modification<DataType,T>> modificationList;

		ModificationManager(multilevelcontext::MultiLevelContextInformation<DataType> multiLevelContext_,
												cosmology::CosmologicalParameters<T> cosmology_,
												fields::OutputField<DataType>* outputField_):
				outputField(outputField_), underlying(multiLevelContext_), cosmology(cosmology_){

		}

		auto calculateValuebyName(std::string name_){

		}

		void addConstrainToList(){

		}

		void applyModifications(){

		}

	private:
		Modification<DataType,T> parseModificationName(std::string name_){
			if ((strcasecmp(name_.c_str(), "overdensity") == 0)){
				auto modif = OverdensityModification();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "potential") == 0)) {
				auto modif = PotentialModification();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "lx") == 0)) {
				auto modif = LxModification();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "ly") == 0)) {
				auto modif = LyModification();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "lz") == 0)) {
				auto modif = LzModification();
				return modif;
			} else{
				std::runtime_error(name_ + "" + "is an unknown modification name");
			}
		}

		std::vector<LinearModification<DataType,T>> createLinearList(){

		}

//		std::vector<QuadraticModification<DataType,T>> createQuadraticList(){
//		}

		void applyLinearModif(){

		}

		void applyQuadModif(){

		}




	};




}






#endif
