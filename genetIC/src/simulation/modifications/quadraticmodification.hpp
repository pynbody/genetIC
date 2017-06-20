#ifndef IC_QUADRATICMODIFICATION_HPP
#define IC_QUADRATICMODIFICATION_HPP

#include <src/simulation/modifications/modification.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class QuadraticModification : public Modification<DataType, T> {
	private:
		int initNumberSteps;
		T targetPrecision;

	public:
		QuadraticModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
													const cosmology::CosmologicalParameters<T> &cosmology_) :
				Modification<DataType, T>(underlying_, cosmology_) {
			this->order = 2;
		};

		int getInitNumberSteps() {
			return this->initNumberSteps;
		}

		T getTargetPrecision(){
			return this->targetPrecision;
		}

		void setInitNumberSteps(int initNumberSteps_){
			this->initNumberSteps = initNumberSteps_;
		}

		void setTargetPrecision(T targetPrecision_){
			this->targetPrecision = targetPrecision_;
		}


		T calculateCurrentValue(fields::MultiLevelField<DataType>*  field) override {

			if(!field->isFourierOnAllLevels()){
				field->toFourier();
			}

			auto pushedField = pushMultiLevelFieldThroughMatrix(*field);
			pushedField.toFourier();

			T value = pushedField.euclidianInnerProduct(*field).real();
			return value;
		}

		fields::MultiLevelField<DataType> pushMultiLevelFieldThroughMatrix(const fields::MultiLevelField<DataType> &field ){

			std::vector<std::shared_ptr<fields::Field<DataType, T>>> oneLevelFieldsVector;
			for (size_t level = 0; level < this->underlying.getNumLevels(); ++level){
				auto pushedOneLevel = pushOneLevelFieldThroughMatrix(field.getFieldForLevel(level));
				oneLevelFieldsVector.push_back(pushedOneLevel);
			}

			fields::MultiLevelField<DataType> multiLevelPushed = fields::MultiLevelField<DataType>(this->underlying,
																																														 oneLevelFieldsVector);
			return multiLevelPushed;
		}

	protected:
		virtual std::shared_ptr<fields::Field<DataType, T>> pushOneLevelFieldThroughMatrix(const fields::Field<DataType, T> &/* field */) = 0;
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class FilteredVarianceModification : public QuadraticModification<DataType, T> {

	public:

		FilteredVarianceModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
		const cosmology::CosmologicalParameters<T> &cosmology_) :
		QuadraticModification<DataType, T>(underlying_, cosmology_) {}

		void setFilterScale(T scale_){
			this->scale =scale_;
		}

		std::shared_ptr<fields::Field<DataType, T>> pushOneLevelFieldThroughMatrix(const fields::Field<DataType, T> &field) override {

			fields::Field<DataType, T> pushedField = fields::Field<DataType, T>(field);

			pushedField.toReal();
			windowOperator(pushedField);

			pushedField.toFourier();
			filterOperator(pushedField);

			pushedField.toReal();
			varianceOperator(pushedField);

			pushedField.toFourier();
			filterOperator(pushedField);

			pushedField.toReal();
			windowOperator(pushedField);

			return std::make_shared<fields::Field<DataType, T>>(pushedField);
		}

	private:
		T scale;

		void windowOperator(fields::Field<DataType, T> &field) {

			assert(! field.isFourier()); // Windowing is done in real space

			std::vector<DataType> &fieldData = field.getDataVector();

			for (size_t i = 0; i < fieldData.size(); ++i) {
				// If cell is not a flagged cell, zero it
				if (!(std::binary_search(this->flaggedCells.begin(), this->flaggedCells.end(), i)) ){
					fieldData[i] = 0;}
			}
		}

		void filterOperator(fields::Field<DataType, T> &field){

			assert(field.isFourier());	//Filtering must be done in Fourier space

			T k_cut = 2 * M_PI / scale;

			auto highPassFermi = filters::ComplementaryFilterAdaptor<filters::LowPassFermiFilter<T>>(k_cut);

			field.applyFilter(highPassFermi);
		}

		void varianceOperator(fields::Field<DataType, T> &field){

			assert(! field.isFourier()); // Variance is calculated in real space

			windowOperator(field);

			std::vector<DataType> &fieldData = field.getDataVector();
			size_t regionSize = this->flaggedCells.size();

			// Calculate mean value in flagged region
			T sum = 0;
			for (size_t i = 0; i < regionSize; i++) {
				 sum += fieldData[this->flaggedCells[i]];
			}

			for (size_t i = 0; i < regionSize; i++) {
				fieldData[this->flaggedCells[i]] -= sum / regionSize;
				fieldData[this->flaggedCells[i]] /= regionSize;
			}
		}

	};
}



#endif
