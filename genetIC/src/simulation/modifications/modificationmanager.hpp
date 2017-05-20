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
		const cosmology::CosmologicalParameters<T> &cosmology;
		std::vector<LinearModification<DataType,T>*> modificationList;

		ModificationManager(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext_,
												cosmology::CosmologicalParameters<T> cosmology_,
												fields::OutputField<DataType>* outputField_):
				outputField(outputField_), underlying(multiLevelContext_), cosmology(cosmology_){

		}

		T calculateCurrentValueByName(std::string name_){

			LinearModification<DataType,T>* modification = getModificationFromName(name_);
			T value = modification->calculateCurrentValue(outputField, underlying);
			//TODO Make sure modification is freed in memory
			return value;
		}

		void addModificationToList(std::string name_, std::string type_ , T target_){
			//TODO Either allow multiple arguments or create a different function for lin and quad
			bool relative = isRelative(type_);
			LinearModification<DataType,T>* modification = getModificationFromName(name_);
			T value = modification->calculateCurrentValue(outputField, underlying);

			T target = target_;
			if(relative) target *= value;

			modification->setTarget(target);

			modificationList.push_back(std::move(modification));
		}

		void applyModifications(){
			applyLinearModif();

		}

	private:
		LinearModification<DataType,T>* getModificationFromName(std::string name_){
			if ((strcasecmp(name_.c_str(), "overdensity") == 0)){
				return new OverdensityModification<DataType,T>(cosmology);
			} else if ((strcasecmp(name_.c_str(), "potential") == 0)) {
				return new PotentialModification<DataType,T>(cosmology);
			} else if ((strcasecmp(name_.c_str(), "lx") == 0)) {
				return new AngMomentumModification<DataType,T>(cosmology, 0);
			} else if ((strcasecmp(name_.c_str(), "ly") == 0)) {
				return new AngMomentumModification<DataType,T>(cosmology, 1);
			} else if ((strcasecmp(name_.c_str(), "lz") == 0)) {
				return new AngMomentumModification<DataType,T>(cosmology, 2);
			} else{
				std::runtime_error(name_ + "" + "is an unknown modification name");
				return NULL;
			}
		}

//		std::vector<LinearModification<DataType,T>> createLinearList(){
//// TODO Implement a way to extract linear modif in a list.
//			return (std::vector<LinearModification<DataType,T>*>) modificationList;
//		}

		void applyLinearModif(){

			std::vector<fields::ConstraintField<DataType>> alphas;
			std::vector<T> targets, existing_values;

			std::cout << "Extract in vectors" << std::endl;
			for (size_t i = 0; i < modificationList.size(); i++) {

				std::cout << "Extract covectors" << std::endl;
				alphas.push_back(std::move(modificationList[i]->calculateCovectorOnAllLevels(underlying)));
				std::cout << "Extract target " << std::endl;
				targets.push_back(modificationList[i]->getTarget());
				std::cout << "Extract values" << std::endl;
				existing_values.push_back(modificationList[i]->calculateCurrentValue(outputField,underlying));

				std::cout << "target = " << targets[i] << std::endl;
				std::cout << "current = " << existing_values[i] << std::endl;
			}


			std::cout << "Starting ortho" << std::endl;
			orthonormaliseModifications(alphas, targets, existing_values);
			std::cout << "Finishing ortho" << std::endl;
			outputField->toFourier();

			std::cout << "Applying constraints" << std::endl;
			for (size_t i = 0; i < alphas.size(); i++) {
				fields::MultiLevelField<DataType> &alpha_i = alphas[i];
				std::cout << "target = " << targets[i] << std::endl;
				std::cout << "current = " << existing_values[i] << std::endl;
				auto dval_i = targets[i] - existing_values[i];

				// Constraints are essentially covectors. Convert to a vector, with correct weighting/filtering.
				// See discussion of why it's convenient to regard constraints as covectors in multilevelfield.hpp
				alpha_i.convertToVector();

				alpha_i.toFourier(); // almost certainly already is in Fourier space, but just to be safe
				outputField->addScaled(alpha_i, dval_i);
			}
		}

		void applyQuadModif(){
			// TODO Implement quadratic algorithm
			/*
			 * for i<N_steps
			 * add linearised quadratic to alphas
			 * orthonormalise constraints
			 * apply linear constraints
			 * repeat
			 */
		}
		void orthonormaliseModifications(std::vector<fields::ConstraintField<DataType>> &alphas,
																		 std::vector<T> &targets, std::vector<T> &existing_values) {
			/* Constraints need to be orthonormal before applying (or alternatively one would need an extra matrix
			 * manipulation on them, as in the original HR91 paper, which boils down to the same thing). */

			using namespace tools::numerics;

			size_t n = alphas.size();
			std::cout << n << std::endl;
			size_t done = 0;

			size_t nCells = underlying.getNumCells();

			// It can be helpful to see a summary of the constraints being applied, with the starting and target values
			std::cout << "v0=" << real(existing_values) << std::endl;
			std::cout << "v1=" << real(targets) << std::endl;

			// Store transformation matrix, just for display purposes
			std::vector<std::vector<std::complex<T>>> t_matrix(n, std::vector<std::complex<T>>(n, 0));

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
					pb.setProgress(((float) done * 2) / (n * (1 + n)));
					T result = alpha_i.innerProduct(alpha_j).real();
					std::cout << "result=" << result << std::endl;

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
				T norm = sqrt(alpha_i.innerProduct(alpha_i).real());
				std::cout << "norm=" << norm << std::endl;
				std::cout << "End Normalize" << std::endl;

				alpha_i /= norm;
				targets[i] /= norm;

				std::cout << "Normalise matrix" << std::endl;
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

		bool isRelative(std::string type){
			bool relative = false;
			if (strcasecmp(type.c_str(), "relative") == 0) {
				relative = true;
			} else if (strcasecmp(type.c_str(), "absolute") != 0) {
				throw std::runtime_error("Modification type must be either 'relative' or 'absolute'");
			}
			return relative;
		}

		bool isLinear(std::string name_){
			bool linear = false;
			if (strcasecmp(name_.c_str(), "overdensity") == 0 || strcasecmp(name_.c_str(), "potential") == 0 ||
					strcasecmp(name_.c_str(), "lx") == 0 || strcasecmp(name_.c_str(), "ly") == 0
					|| strcasecmp(name_.c_str(), "lz") == 0 ){
				linear = true;
			} else{
				throw std::runtime_error(name_ + "" + " is not an implemented linear modification");
			}
			return linear;
		}

		bool isQuadratic(std::string name_){
			bool linear = false;
			if (strcasecmp(name_.c_str(), "variance") == 0){
				linear = true;
			} else{
				throw std::runtime_error(name_ + "" + " is not an implemented quadratic modification");
			}
			return linear;
		}

		bool areInputCorrect(){
			//TODO Method to check variable arguments are correct in addModiftolist.
			return false;
		}


	};




}






#endif
