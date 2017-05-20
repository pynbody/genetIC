#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/field/multilevelfield.hpp>
#include <src/simulation/modifications/modification.hpp>
#include <src/cosmology/parameters.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearModification : public Modification<DataType, T> {

	public:

		LinearModification(const cosmology::CosmologicalParameters<T> &cosmology_):Modification<DataType,T>(cosmology_){};

		DataType calculateCurrentValue(fields::MultiLevelField<DataType>* field,
														multilevelcontext::MultiLevelContextInformation<DataType> &underlying) override {
			//TODO Decide whether context should be extracted from field or passed as an extra argument
			auto covector = calculateCovectorOnAllLevels(underlying);
			DataType val = (field->innerProduct(covector)).real();
			return val;
		}

		fields::ConstraintField<DataType> calculateCovectorOnAllLevels(multilevelcontext::MultiLevelContextInformation<DataType> &underlying) {

			int level = underlying.getNumLevels() - 1;

			using tools::numerics::operator/=;
			auto highResConstraint = calculateCovectorOnOneLevel(this->cosmology, underlying.getGridForLevel(level));

			if (level != 0)
				highResConstraint.getDataVector() /= underlying.getWeightForLevel(level);

			auto covector = underlying.generateMultilevelFromHighResField(std::move(highResConstraint));
			return covector;
		}


	protected:
		virtual fields::Field<DataType, T> calculateCovectorOnOneLevel(const cosmology::CosmologicalParameters<T> &cosmology,
																																	 grids::Grid<DataType> &grid) = 0;
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class OverdensityModification : public LinearModification<DataType, T> {
	public:

		OverdensityModification(const cosmology::CosmologicalParameters<T> &cosmology_):LinearModification<DataType,T>(cosmology_){};

		fields::Field<DataType, T> calculateCovectorOnOneLevel(const cosmology::CosmologicalParameters<T> & /*&cosmology*/,
																													 grids::Grid<DataType> &grid) override {

			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			std::vector<DataType> outputData = outputField.getDataVector();
			std::vector<size_t> particleArray;
			grid.getFlaggedCells(particleArray);

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
		PotentialModification(const cosmology::CosmologicalParameters<T> &cosmology_):LinearModification<DataType,T>(cosmology_){};

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
		AngMomentumModification(const cosmology::CosmologicalParameters<T> &cosmology_, int direction_):LinearModification<DataType,T>(cosmology_){
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

