#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/field/multilevelfield.hpp>
#include <src/simulation/modifications/modification.hpp>
#include <src/cosmology/parameters.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearModification : public Modification<DataType, T> {
	protected:
		fields::ConstraintField<DataType> covector;

		bool isCovectorInitialised(){

		}

	public:
		T calculateCurrentValue(fields::MultiLevelField<DataType> field,
														cosmology::CosmologicalParameters<T> &cosmology_) override {
			calculateCovectorOnAllLevels(field.getContext(), cosmology_);
			T val = covector.innerProduct(field);
			return val;
		}

		void calculateCovectorOnAllLevels(multilevelcontext::MultiLevelContextInformation<DataType> &underlying,
																			cosmology::CosmologicalParameters<T> &cosmology) {

			int level = underlying.getNumLevels() - 1;

			using tools::numerics::operator/=;
			auto highResConstraint = calculateCovectorOnOneLevel(cosmology, underlying.getGridForLevel(level));

			if (level != 0)
				highResConstraint.getDataVector() /= underlying.getWeightForLevel(level);

			covector = std::move(underlying.generateMultilevelFromHighResField(std::move(highResConstraint)));
		}




	protected:

		virtual fields::Field<DataType, T> calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology,
																																	 grids::Grid<DataType> &grid) = 0;

	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class OverdensityModification : public LinearModification<DataType, T> {
	public:
		fields::Field<DataType, T> calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology,
																													 grids::Grid<DataType> &grid) override {

			fields::Field<DataType, T> outputField = fields::Field(grid,false);
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
		void
		calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology, grids::Grid<DataType> &grid) override {

		}
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LxModification : public LinearModification<DataType, T> {
	public:
		void
		calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology, grids::Grid<DataType> &grid) override {

		}
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LyModification : public LinearModification<DataType, T> {
	public:
		void
		calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology, grids::Grid<DataType> &grid) override {

		}
	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LzModification : public LinearModification<DataType, T> {
	public:
		void
		calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology, grids::Grid<DataType> &grid) override {

		}

	};

}

#endif

