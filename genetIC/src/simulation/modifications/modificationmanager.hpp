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
			auto modification = getModificationFromName(name_);
			auto value = modification.calculateCurrentValue(*outputField);
			return value;
		}

		void addModificationToList(){

		}

		void applyModifications(){
			applyLinearModif();

		}

	private:
		Modification<DataType,T> getModificationFromName(std::string name_){
			if ((strcasecmp(name_.c_str(), "overdensity") == 0)){
				auto modif = OverdensityModification<DataType,T>();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "potential") == 0)) {
				auto modif = PotentialModification<DataType,T>();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "lx") == 0)) {
				auto modif = LxModification<DataType,T>();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "ly") == 0)) {
				auto modif = LyModification<DataType,T>();
				return modif;
			} else if ((strcasecmp(name_.c_str(), "lz") == 0)) {
				auto modif = LzModification<DataType,T>();
				return modif;
			} else{
				std::runtime_error(name_ + "" + "is an unknown modification name");
			}
		}

		std::vector<LinearModification<DataType,T>> createLinearList(){

		}

//		std::vector<QuadraticModification<DataType,T>> createQuadraticList(){
//		}

		void orthonormaliseModifications(std::vector<fields::ConstraintField<DataType>> alphas,
																			std::vector<T> targets, std::vector<T> existing_values) {
			/* Constraints need to be orthonormal before applying (or alternatively one would need an extra matrix
			 * manipulation on them, as in the original HR91 paper, which boils down to the same thing). */

			using namespace tools::numerics;

			size_t n = alphas.size();
			size_t done = 0;

			size_t nCells = underlying.getNumCells();

			// It can be helpful to see a summary of the constraints being applied, with the starting and target values
			std::cout << "v0=" << real(existing_values) << std::endl;
			std::cout << "v1=" << real(targets) << std::endl;

			// Store transformation matrix, just for display purposes
			std::cout << "Define matrix" << std::endl;
			std::vector<std::vector<std::complex<T>>> t_matrix(n, std::vector<std::complex<T>>(n, 0));

			std::cout << "First loop" <<  std::endl;
			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < n; j++) {
					if (i == j) t_matrix[i][j] = 1.0; else t_matrix[i][j] = 0.0;
				}
			}

			// Gram-Schmidt orthogonalization in-place

			tools::progress::ProgressBar pb("orthogonalizing constraints");

			for (size_t i = 0; i < n; i++) {
				auto &alpha_i = alphas[i];
				for (size_t j = 0; j < i; j++) {
					auto &alpha_j = alphas[j];
					std::cout << alpha_j.isFourierOnAllLevels() << std::endl;
					pb.setProgress(((float) done * 2) / (n * (1 + n)));
					T result = alpha_i.innerProduct(alpha_j).real();

					alpha_i.addScaled(alpha_j, -result);


					// update display matrix accordingly such that at each step
					// our current values array is equal to t_matrix . input_values
					for (size_t k = 0; k <= j; k++)
						t_matrix[i][k] -= result * t_matrix[j][k];

					// update constraining value
					targets[i] -= result * targets[j];
					done += 1; // one op for each of the orthognalizing constraints
				}

				// normalize
				std::cout << "Normalize" << std::endl;
				std::cout << alpha_i.isFourierOnAllLevels() << std::endl;
				T norm = sqrt(alpha_i.innerProduct(alpha_i).real());
				std::cout << "End Normalize" << std::endl;

				alpha_i /= norm;
				targets[i] /= norm;

				std::cout << "Final normalization loop" << std::endl;
				for (size_t j = 0; j < n; j++) {
					t_matrix[i][j] /= norm;
				}
				done += 1; // one op for the normalization

			}

			// The existing values are now invalid because the constraint vectors have been changed.
			// Re-calculate them.
			// TODO: rather than recalculate them, it would save CPU time to update them alongside the vectors/target values.
			// This is actually pretty trivial but needs to be carefully tested.
			std::cout << "Recalculate existing" << std::endl;
			existing_values.clear();
			for (size_t i = 0; i < n; i++) {
				auto &alpha_i = alphas[i];
				existing_values.push_back(alpha_i.innerProduct(*outputField).real());
			}


			// Finally display t_matrix^{dagger} t_matrix, which is the matrix allowing one
			// to form chi^2: Delta chi^2 = d_1^dagger t_matrix^dagger t_matrix d_1 -
			// d_0^dagger t_matrix^dagger t_matrix d_0.

			std::cout << std::endl << "chi2_matr = [";
			for (size_t i = 0; i < n; i++) {
				std::cout << "[";
				for (size_t j = 0; j < n; j++) {
					std::complex<T> r = 0;
					for (size_t k = 0; k < n; k++) {
						r += std::conj(t_matrix[k][i]) * t_matrix[k][j];
					}
					std::cout << std::real(r) / nCells;
					if (j < n - 1) std::cout << ",";
				}
				std::cout << "]";
				if (i < n - 1) std::cout << "," << std::endl;
			}
			std::cout << "]" << std::endl;
		}

		void applyLinearModif(){
			auto linearList = createLinearList();

			// TODO Translate linear list into vectors of alphas, targets and existing
			std::vector<fields::ConstraintField<DataType>> alphas = createAlphaList();
			std::vector<T> targets = createTargetsList();
			std::vector<T> existing_values = createExistingList();


			std::cout << "Starting ortho" << std::endl;
			orthonormaliseModifications(alphas, targets, existing_values);
			std::cout << "Finishing ortho" << std::endl;
			std::cout << "Pointer to Fourier" << std::endl;
			outputField->toFourier();

			std::cout << "Applying constraints" << std::endl;
			for (size_t i = 0; i < alphas.size(); i++) {
				fields::MultiLevelField<DataType> &alpha_i = alphas[i];
				auto dval_i = targets[i] - existing_values[i];

				// Constraints are essentially covectors. Convert to a vector, with correct weighting/filtering.
				// See discussion of why it's convenient to regard constraints as covectors in multilevelfield.hpp
				alpha_i.convertToVector();

				alpha_i.toFourier(); // almost certainly already is in Fourier space, but just to be safe
				outputField->addScaled(alpha_i, dval_i);
			}
		}

		void applyQuadModif(){
			/*
			 * for i<N_steps
			 * add linearised quadratic to alphas
			 * orthonormalise constraints
			 * apply linear constraints
			 * repeat
			 */
		}




	};




}






#endif
