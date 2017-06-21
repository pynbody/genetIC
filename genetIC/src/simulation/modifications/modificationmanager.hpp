#ifndef IC_MODIFICATIONMANAGER_HPP
#define IC_MODIFICATIONMANAGER_HPP

#include <src/simulation/modifications/linearmodification.hpp>
#include <src/simulation/modifications/quadraticmodification.hpp>

//! Deals with the creation of genetically modified fields
namespace modifications{

	//! Keep track of modifications to be applied. Construct the modified field from these modifications.
	template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
	class ModificationManager {

	public:
		fields::OutputField<DataType>* outputField;																			/*!< Will become the modified field */
		multilevelcontext::MultiLevelContextInformation<DataType> &underlying;					/*!< Grid context in which modifications take place */
		const cosmology::CosmologicalParameters<T> &cosmology;													/*!< Cosmology context in which modifications take place */

		std::vector<std::shared_ptr<LinearModification<DataType,T>>> linearModificationList;  /*!< Modifications to be applied */
		std::vector<std::shared_ptr<QuadraticModification<DataType,T>>> quadraticModificationList;

		ModificationManager(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext_,
												const cosmology::CosmologicalParameters<T> &cosmology_,
												fields::OutputField<DataType>* outputField_):
				outputField(outputField_), underlying(multiLevelContext_), cosmology(cosmology_){}

		//! Calculate existing value of the quantity defined by name
		T calculateCurrentValueByName(std::string name_){

			std::shared_ptr<Modification<DataType,T>> modification = getModificationFromName(name_);
			T value = modification->calculateCurrentValue(outputField);
			return value;
		}

		T calculateVariance(T scale_){

			std::shared_ptr<FilteredVarianceModification<DataType,T>> modification =
					make_shared<FilteredVarianceModification<DataType,T>>(underlying, cosmology);
			modification->setFilterScale(scale_);
			T value = modification->calculateCurrentValue(outputField);
			return value;
		}

		/*!
				\param type_ Modification can be relative to existing value or absolute
				\param target_ Absolute target or factor by which the existing will be multiplied
			*/
		void addModificationToList(std::string name_, std::string type_ , T target_){
      
			std::shared_ptr<Modification<DataType,T>> modification = getModificationFromName(name_);

			bool relative = isRelative(type_);
			T target = target_;

			if(relative){
				T value = modification->calculateCurrentValue(outputField);
				target *= value;
			}

			modification->setTarget(target);

			if(modification->getOrder() == 1 ){
				linearModificationList.push_back(std::dynamic_pointer_cast<LinearModification<DataType,T>>(modification));
			} else if (modification->getOrder() == 2){
				quadraticModificationList.push_back(std::dynamic_pointer_cast<QuadraticModification<DataType,T>>(modification));
			} else{
				throw std::runtime_error( " Could not add modification to list");
			}
		}


		void addQuadModificationToList(std::string /*name_*/, std::string type_ , T target_, int initNsteps_, T precision_, T scale_ = 0){

			std::shared_ptr<FilteredVarianceModification<DataType,T>> modification =
									make_shared<FilteredVarianceModification<DataType,T>>(underlying, cosmology);

			modification->setInitNumberSteps(initNsteps_);
			modification->setTargetPrecision(precision_);
			modification->setFilterScale(scale_);

			bool relative = isRelative(type_);
			T target = target_;

			if(relative){
				T value = modification->calculateCurrentValue(outputField);
				target *= value;
			}

			modification->setTarget(target);

			quadraticModificationList.push_back(modification);
		}

		//! Construct the modified field with all modifications present in the modification list
		void applyModifications(){
			applyLinearModif();

		}

		//! Calculate and displays the covariance matrix of the modification covectors
		/*!
		 * Mostly for debugging purposes
 		*/
		void print_covariance() {

			tools::progress::ProgressBar pb("calculating covariance");

			std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas;

			for (size_t i = 0; i < linearModificationList.size(); i++) {
				alphas.push_back(linearModificationList[i]->getCovector());
			}

			size_t n = alphas.size();
			size_t done = 0;
			std::vector<std::vector<std::complex<T>>> c_matrix(n, std::vector<std::complex<T>>(n, 0));

			for (size_t i = 0; i < n; i++) {
				auto &alpha_i = *(alphas[i]);
				for (size_t j = 0; j <= i; j++) {
					pb.setProgress(((float) done * 2) / (n * (1 + n)));
					auto &alpha_j = *(alphas[j]);
					c_matrix[i][j] = alpha_i.innerProduct(alpha_j).real();
					c_matrix[j][i] = c_matrix[i][j];
					done += 1;
				}
			}

			std::cout << std::endl << "cov_matr = [";
			for (size_t i = 0; i < n; i++) {
				std::cout << "[";
				for (size_t j = 0; j < n; j++) {
					std::cout << std::real(c_matrix[i][j]);
					if (j < n - 1) std::cout << ",";
				}
				std::cout << "]";
				if (i < n - 1) std::cout << "," << std::endl;
			}
			std::cout << "]" << std::endl;
		}

		void clearModifications(){
			std::cout << "Clearing modification list" << std::endl;
			linearModificationList.clear();
		}



	private:
		std::shared_ptr<Modification<DataType,T>> getModificationFromName(std::string name_){
			try{
				auto modification = getLinearModificationFromName(name_);
				return modification;
			} catch( UnknownModificationException &e) {
				// If modification is unknwon, it might be quadratic so swallow exception for now.
				try {
					auto modification = getQuadraticModificationFromName(name_);
					return modification;
				} catch (UnknownModificationException &e2) {
					throw e2;
				}
			}

		}


		std::shared_ptr<LinearModification<DataType,T>> getLinearModificationFromName(std::string name_){
			if ((strcasecmp(name_.c_str(), "overdensity") == 0)){
				return make_shared<OverdensityModification<DataType,T>>(underlying, cosmology);
			} else if ((strcasecmp(name_.c_str(), "potential") == 0)) {
				return make_shared<PotentialModification<DataType,T>>(underlying, cosmology);
			} else if ((strcasecmp(name_.c_str(), "lx") == 0)) {
				return make_shared<AngMomentumModification<DataType,T>>(underlying, cosmology, 0);
			} else if ((strcasecmp(name_.c_str(), "ly") == 0)) {
				return make_shared<AngMomentumModification<DataType,T>>(underlying, cosmology, 1);
			} else if ((strcasecmp(name_.c_str(), "lz") == 0)) {
				return make_shared<AngMomentumModification<DataType,T>>(underlying, cosmology, 2);
			} else{
				throw UnknownModificationException(name_ + " " + "is an unknown modification name");
			}
		}

		std::shared_ptr<QuadraticModification<DataType,T>> getQuadraticModificationFromName(std::string name_){
			if ((strcasecmp(name_.c_str(), "variance") == 0)){
				return make_shared<FilteredVarianceModification<DataType,T>>(underlying, cosmology);
			}  else{
				throw UnknownModificationException(name_ + " " + "is an unknown modification name");
			}
		}


		/*!
		 * Linear modifications are applied by orthonormalisation and adding the
		 * correction term. See Roth et al 2016 for details
		 */
		void applyLinearModif(){

			std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas;
			std::vector<T> targets, existing_values;

			// Extract A, b and A*delta_0 from modification list
			for (size_t i = 0; i < linearModificationList.size(); i++) {
				alphas.push_back(linearModificationList[i]->getCovector());
				targets.push_back(linearModificationList[i]->getTarget());
				existing_values.push_back(linearModificationList[i]->calculateCurrentValue(outputField));
			}


			orthonormaliseModifications(alphas, targets, existing_values);
			outputField->toFourier();

			for (size_t i = 0; i < alphas.size(); i++) {
				auto &alpha_i = *(alphas[i]);
				auto dval_i = targets[i] - existing_values[i];

				// Constraints are essentially covectors. Convert to a vector, with correct weighting/filtering.
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

			size_t numberQuadraticModifs = quadraticModificationList.size();
			//TODO Check quadratic are indepednent ?

			for (size_t i=0; i<numberQuadraticModifs; i++){
				std::shared_ptr<QuadraticModification<DataType,T>> modif_i = quadraticModificationList[i];

				T starting_value = modif_i->calculateCurrentValue(this->outputField);
				T overall_target = modif_i->getTarget();
				int n_steps = modif_i->getInitNumberSteps();

				std::vector<T> quad_targets = tools::linspace(starting_value, overall_target, n_steps);

				performIteration(modif_i, quad_targets, n_steps);

				n_steps = calculateCorrectNumberSteps(modif_i, n_steps);
				quad_targets = tools::linspace(starting_value, overall_target, n_steps);

				performIteration(modif_i, quad_targets, n_steps);

			}
		}

		void performIteration(std::shared_ptr<QuadraticModification<DataType,T>> modif_, std::vector<T> quad_targets_, int n_steps_){

			std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas;
			std::vector<T> linear_targets, linear_existing_values;

			// Extract A, b and A*delta_0 from modification list
			for (size_t i = 0; i < linearModificationList.size(); i++) {
				alphas.push_back(linearModificationList[i]->getCovector());
				linear_targets.push_back(linearModificationList[i]->getTarget());
				linear_existing_values.push_back(linearModificationList[i]->calculateCurrentValue(outputField));
			}


			for (size_t l=0; l<n_steps_;l++){


				auto pushedField = modif_->pushMultiLevelFieldThroughMatrix(*outputField);

				alphas.push_back(pushedField); //Add pushed field to linear covectors
				linear_targets.push_back(0);
				linear_existing_values.push_back(0);

				orthonormaliseModifications(alphas, linear_targets, linear_existing_values);

				alphas.pop_back();
				linear_targets.pop_back();
				linear_existing_values.pop_back();

				applyLinearModif();

				T quadratic_existing = modif_->calculateCurrentValue(outputField);
				T multiplier = (1/2) * (quad_targets_[l] - quadratic_existing);

				pushedField->convertToVector();
				outputField->addScaled(*pushedField, multiplier);
			}
		}

		int calculateCorrectNumberSteps(std::shared_ptr<QuadraticModification<DataType,T>> modif_, int previous_n_steps){
			T targeted_precision = modif_->getTargetPrecision();
			T current_value = modif_->calculateCurrentValue(outputField);
			T overall_target = modif_->getTarget();

			T residual = overall_target - current_value;
			T achieved_precision = residual / overall_target;

			if(achieved_precision < targeted_precision){
				return 0;
			} else {
				T scaling = previous_n_steps * sqrt(achieved_precision/targeted_precision);
				return (int) std::ceil(scaling);
			}

		}

		//! Graam-Schmidt procedure to orthonormalise the modification covectors
		void orthonormaliseModifications(std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas,
																		 std::vector<T> &targets, std::vector<T> &existing_values) {

			using namespace tools::numerics;

			size_t n = alphas.size();
			size_t done = 0;

			size_t nCells = underlying.getNumCells();

			// Summary of modifs to be applied
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
			tools::progress::ProgressBar pb("orthogonalizing modifications");

			for (size_t i = 0; i < n; i++) {
				auto &alpha_i = *(alphas[i]);
				for (size_t j = 0; j < i; j++) {
					auto &alpha_j = *(alphas[j]);
					pb.setProgress(((float) done * 2) / (n * (1 + n)));
					T result = alpha_i.innerProduct(alpha_j).real();

					alpha_i.addScaled(alpha_j, -result);


					// update display matrix accordingly such that at each step
					// our current values array is equal to t_matrix.
					for (size_t k = 0; k <= j; k++)
						t_matrix[i][k] -= result * t_matrix[j][k];

					// update constraining value
					targets[i] -= result * targets[j];
					done += 1; // one op for each of the orthognalizing modifs
				}

				// normalize
				T norm = sqrt(alpha_i.innerProduct(alpha_i).real());

				alpha_i /= norm;
				targets[i] /= norm;

				for (size_t j = 0; j < n; j++) {
					t_matrix[i][j] /= norm;
				}
				done += 1; // one op for the normalization

			}

			// The existing values are now invalid because the covectors have been changed.
			// TODO: rather than recalculate them, it would save CPU time to update them alongside the vectors/target values.
			existing_values.clear();
			for (size_t i = 0; i < n; i++) {
				auto &alpha_i = *(alphas[i]);
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

		bool areInputCorrect(){
			//TODO Method to check variable arguments are correct in addModiftolist.
			return false;
		}


	};




}





#endif
