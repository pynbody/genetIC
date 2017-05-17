#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/field/multilevelfield.hpp>
#include <src/simulation/modifications/modification.hpp>
#include <src/cosmology/parameters.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearModification : public Modification<DataType, T> {
	private:
		fields::ConstraintField<DataType> covector;

	public:
		T calculateCurrentValue(fields::MultiLevelField<DataType> field) override {
			T val = covector.innerProduct(field);
			return val;
		}

		void calculateCovectorOnAllLevels(multilevelcontext::MultiLevelContextInformation<DataType> &underlying,
																			cosmology::CosmologicalParameters<T> &cosmology) {

		}

		virtual void
		calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology, grids::Grid<DataType> &grid) = 0;

	};

	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class OverdensityModification : public LinearModification<DataType, T> {
	public:
		void
		calculateCovectorOnOneLevel(cosmology::CosmologicalParameters<T> &cosmology, grids::Grid<DataType> &grid) override {

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

