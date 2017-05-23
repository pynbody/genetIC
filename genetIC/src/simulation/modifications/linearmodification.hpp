#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/modifications/modification.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearModification : public Modification<DataType, T> {

	public:

		//TODO underlying could be a const reference depending on what Andrew says
		LinearModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
											 const cosmology::CosmologicalParameters<T> &cosmology_): Modification<DataType,T>(underlying_, cosmology_){};

		DataType calculateCurrentValue(fields::MultiLevelField<DataType>* field) override {
			//TODO Decide whether context should be extracted from field or passed as an extra argument
			auto covector = calculateCovectorOnAllLevels();
			covector.toFourier();
			DataType val = covector.innerProduct(*field).real();
			return val;
		}

		fields::ConstraintField<DataType> calculateCovectorOnAllLevels() {

			size_t level = this->underlying.getNumLevels() - 1;

			using tools::numerics::operator/=;
			auto highResConstraint = calculateCovectorOnOneLevel(this->cosmology, this->underlying.getGridForLevel(level));

			if (level != 0) {
				highResConstraint.getDataVector() /= this->underlying.getWeightForLevel(level);
			}

			auto covector = this->underlying.generateMultilevelFromHighResField(std::move(highResConstraint));
			covector.toFourier();
			return std::move(covector);
		}


	protected:
		virtual fields::Field<DataType, T> calculateCovectorOnOneLevel(const cosmology::CosmologicalParameters<T> &cosmology,
																																	 grids::Grid<DataType> &grid) = 0;
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class OverdensityModification : public LinearModification<DataType, T> {
	public:

		OverdensityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
														const cosmology::CosmologicalParameters<T> &cosmology_): LinearModification<DataType,T>(underlying_, cosmology_){};

		fields::Field<DataType, T> calculateCovectorOnOneLevel(const cosmology::CosmologicalParameters<T> & /*&cosmology*/,
																													 grids::Grid<DataType> &grid) override {

			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			std::vector<DataType> &outputData = outputField.getDataVector();
			std::vector<size_t> particleArray;

			// Get the region of interest and store in particle array
			grid.getFlaggedCells(particleArray);
			for (auto i = particleArray.begin(); i != particleArray.end(); ++i)
				std::cout << *i << ' ';


			T w = 1.0 / particleArray.size();

			for (size_t i = 0; i < grid.size3; ++i) {
				outputData[i] = 0;
			}

			for (size_t i = 0; i < particleArray.size(); i++) {
				outputData[particleArray[i]] += w;
			}

			outputField.toFourier();
			return outputField;
		}
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class PotentialModification : public LinearModification<DataType, T> {
	public:
		PotentialModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
													const cosmology::CosmologicalParameters<T> &cosmology_): LinearModification<DataType,T>(underlying_, cosmology_){};

		fields::Field<DataType, T> calculateCovectorOnOneLevel(const cosmology::CosmologicalParameters<T> & /*&cosmology*/,
																													 grids::Grid<DataType> &grid) override {
			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			return outputField;

		}
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class AngMomentumModification : public LinearModification<DataType, T> {
	public:
		int direction;

	public:
		AngMomentumModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
														const cosmology::CosmologicalParameters<T> &cosmology_, int direction_):
			LinearModification<DataType,T>(underlying_, cosmology_){
			direction = direction_;
		};

		fields::Field<DataType, T> calculateCovectorOnOneLevel(const cosmology::CosmologicalParameters<T> & /*&cosmology*/,
																													 grids::Grid<DataType> &grid) override {
			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			return outputField;
		}
	};

}

#endif

