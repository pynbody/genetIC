#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/modifications/modification.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearModification : public Modification<DataType, T> {

	public:

		LinearModification(const cosmology::CosmologicalParameters<T> &cosmology_):Modification<DataType,T>(cosmology_){};

		DataType calculateCurrentValue(fields::MultiLevelField<DataType>* field,
														multilevelcontext::MultiLevelContextInformation<DataType> &underlying) override {
			//TODO Decide whether context should be extracted from field or passed as an extra argument
			std::cout << "Calculate covector" << std::endl;
			auto covector = calculateCovectorOnAllLevels(underlying);
			std::cout << "Taking inner Product" << std::endl;
			DataType val = (field->innerProduct(covector)).real();
			return val;
		}

		fields::ConstraintField<DataType> calculateCovectorOnAllLevels(multilevelcontext::MultiLevelContextInformation<DataType> &underlying) {

			size_t level = underlying.getNumLevels();

			using tools::numerics::operator/=;
			std::cout << "Calculate highRes" << std::endl;
			auto highResConstraint = calculateCovectorOnOneLevel(this->cosmology, underlying.getGridForLevel(level));
			std::cout << "Leaving OneLevel" << std::endl;

			if (level != 0) {
				std::cout << "divide by weight" << std::endl;
				highResConstraint.getDataVector() /= underlying.getWeightForLevel(level);
			}

			std::cout << "Generate multigrid covector" << std::endl;
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

			std::cout << "Declare one level" << std::endl;
			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			std::cout << "Declare output" << std::endl;
			std::vector<DataType> outputData = outputField.getDataVector();
			std::cout << "Declare partcile array" << std::endl;
			std::vector<size_t> particleArray;
			std::cout << "get Flagged" << std::endl;
			grid.getFlaggedCells(particleArray);

			T w = 1.0 / particleArray.size();

			for (size_t i = 0; i < grid.size3; ++i) {
				outputData[i] = 0;
			}

			for (size_t i = 0; i < particleArray.size(); i++) {
				outputData[particleArray[i]] += w;
			}
			std::cout << "to Fourier" << std::endl;
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

