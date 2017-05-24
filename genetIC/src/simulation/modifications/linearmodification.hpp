#ifndef IC_LINEARMODIFICATION_HPP
#define IC_LINEARMODIFICATION_HPP

#include <src/simulation/modifications/modification.hpp>

namespace modifications {


	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class LinearModification : public Modification<DataType, T> {

	public:

		//TODO underlying could be a const reference depending on what Andrew says
		LinearModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
											 const cosmology::CosmologicalParameters<T> &cosmology_):
				Modification<DataType,T>(underlying_, cosmology_){};

		DataType calculateCurrentValue(fields::MultiLevelField<DataType>* field) override {
			auto covector = calculateCovectorOnAllLevels();
			covector.toFourier();
			DataType val = covector.innerProduct(*field).real();
			return val;
		}

		fields::ConstraintField<DataType> calculateCovectorOnAllLevels() {

			size_t level = this->underlying.getNumLevels() - 1;

			using tools::numerics::operator/=;
			auto highResConstraint = calculateCovectorOnOneLevel(this->underlying.getGridForLevel(level));

			if (level != 0) {
				highResConstraint.getDataVector() /= this->underlying.getWeightForLevel(level);
			}

			auto covector = this->underlying.generateMultilevelFromHighResField(std::move(highResConstraint));
			covector.toFourier();
			return std::move(covector);
		}



	protected:

		virtual fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<DataType> &grid) = 0;




		Coordinate<T> getCentre(grids::Grid<DataType> &grid) {

			T xa, ya, za, xb, yb, zb, x0 = 0., y0 = 0., z0 = 0.;


			std::tie(xa, ya, za) = grid.getCellCentroid(this->flaggedCells[0]);

			for (size_t i = 0; i < this->flaggedCells.size(); i++) {
				std::tie(xb, yb, zb) = grid.getCellCentroid(this->flaggedCells[i]);

				x0 += grid.getWrappedOffset(xb, xa);
				y0 += grid.getWrappedOffset(yb, ya);
				z0 += grid.getWrappedOffset(zb, za);
			}

			x0 /= this->flaggedCells.size();
			y0 /= this->flaggedCells.size();
			z0 /= this->flaggedCells.size();
			x0 += xa;
			y0 += ya;
			z0 += za;

			Coordinate<T> result = Coordinate<T>(x0, y0, z0);
			return result;
		}

		void centralDifference4thOrder(grids::Grid<DataType> &grid, std::vector<DataType> &outputData, long index, int direc, T x0, T y0, T z0) {

			T xp = 0., yp = 0., zp = 0.;
			std::tie(xp, yp, zp) = grid.getCellCentroid(index);

			xp = grid.getWrappedOffset(xp, x0);
			yp = grid.getWrappedOffset(yp, y0);
			zp = grid.getWrappedOffset(zp, z0);

			T c[3] = {0, 0, 0};
			if (direc == 0) {
				c[2] = yp;
				c[1] = -zp;
			} else if (direc == 1) {
				c[0] = zp;
				c[2] = -xp;
			} else if (direc == 2) {
				c[1] = xp;
				c[0] = -yp;
			} else if (direc == 3) {
				T rp = std::sqrt((xp * xp) + (yp * yp) + (zp * zp));
				if (rp != 0) {
					c[0] = xp / rp;
					c[1] = yp / rp;
					c[2] = zp / rp;
				}
			} // radial velocity

			else {
				throw std::runtime_error("Wrong value for parameter 'direc' in function 'cen_deriv4_alpha'");
			}

			for (int di = 0; di < 3; di++) {
				long ind_p1, ind_m1, ind_p2, ind_m2;
				int step1[3] = {0, 0, 0};
				int neg_step1[3] = {0, 0, 0};
				step1[di] = 1;
				neg_step1[di] = -1;

				// N.B. can't wrap - might be on subgrid -> we do it anyway! this will likely break with zoom-in TODO
				//ind_m1=grid.findNextIndNoWrap(index, neg_step1);
				//ind_p1=grid.findNextIndNoWrap(index, step1);
				//ind_m2=grid.findNextIndNoWrap(ind_m1, neg_step1);
				//ind_p2=grid.findNextIndNoWrap(ind_p1, step1);

				ind_m1 = grid.getIndexFromIndexAndStep(index, neg_step1);
				ind_p1 = grid.getIndexFromIndexAndStep(index, step1);
				ind_m2 = grid.getIndexFromIndexAndStep(ind_m1, neg_step1);
				ind_p2 = grid.getIndexFromIndexAndStep(ind_p1, step1);

				T a = -1. / 12. / grid.cellSize, b = 2. / 3. / grid.cellSize;  //the signs here so that L ~ - Nabla Phi

				outputData[ind_m2] += (c[di] * a);
				outputData[ind_m1] += (c[di] * b);
				outputData[ind_p1] += (-c[di] * b);
				outputData[ind_p2] += (-c[di] * a);
			}

		}
	};







	template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
	class OverdensityModification : public LinearModification<DataType, T> {
	public:

		OverdensityModification(multilevelcontext::MultiLevelContextInformation<DataType> &underlying_,
														const cosmology::CosmologicalParameters<T> &cosmology_):
				LinearModification<DataType,T>(underlying_, cosmology_){};

		fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<DataType> &grid) override {

			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			std::vector<DataType> &outputData = outputField.getDataVector();


			T w = 1.0 / this->flaggedCells.size();

			for (size_t i = 0; i < grid.size3; ++i) {
				outputData[i] = 0;
			}

			for (size_t i = 0; i < this->flaggedCells.size(); i++) {
				outputData[this->flaggedCells[i]] += w;
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

		fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<DataType> &grid) override {
			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			std::vector<DataType> &outputData = outputField.getDataVector();

			T w = 1.0 / this->flaggedCells.size();

			for (size_t i = 0; i < grid.size3; ++i) {
				outputData[i] = 0;
			}

			for (size_t i = 0; i < this->flaggedCells.size(); i++) {
				outputData[this->flaggedCells[i]] += w;
			}

			outputField.toFourier();
			densityToPotential(outputField, this->cosmology);
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

			if (direction_ < 0 || direction_ > 2)
				throw std::runtime_error("Angular momentum direction must be 0 (x), 1 (y) or 2 (z)");

			direction = direction_;
		};

		fields::Field<DataType, T> calculateCovectorOnOneLevel(grids::Grid<DataType> &grid) override {



			fields::Field<DataType, T> outputField = fields::Field<DataType, T>(grid,false);
			std::vector<DataType> &outputData = outputField.getDataVector();


			T x0, y0, z0;
			x0 = this->getCentre(grid).x;
			y0 = this->getCentre(grid).y;
			z0 = this->getCentre(grid).z;

			for (size_t i = 0; i < this->flaggedCells.size(); i++) {
				this->centralDifference4thOrder(grid, outputData, this->flaggedCells[i], this->direction, x0, y0, z0);
			}

			outputField.toFourier();

			// The constraint as derived is on the potential. By considering
			// unitarity of FT, we can FT the constraint to get the constraint
			// on the density.
			densityToPotential(outputField, this->cosmology);

			return outputField;
		}
	};

}

#endif

