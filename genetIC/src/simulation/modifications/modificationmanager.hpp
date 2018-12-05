#ifndef IC_MODIFICATIONMANAGER_HPP
#define IC_MODIFICATIONMANAGER_HPP

#include <src/simulation/modifications/linearmodification.hpp>
#include <src/simulation/modifications/quadraticmodification.hpp>

//! Deals with the creation of genetically modified fields
namespace modifications {

  //! Keep track of modifications to be applied. Construct the modified field from these modifications.
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class ModificationManager {

  public:
    std::vector<std::shared_ptr<fields::OutputField<DataType>>> &outputFields;      /*!< All fields to propagate modifications to*/
    multilevelcontext::MultiLevelContextInformation<DataType> &underlying;          /*!< Grid context in which modifications take place */
    const cosmology::CosmologicalParameters<T> &cosmology;                          /*!< Cosmology context in which modifications take place */

    std::vector<std::shared_ptr<LinearModification<DataType, T>>> linearModificationList;  /*!< Modifications to be applied */
    std::vector<std::shared_ptr<QuadraticModification<DataType, T>>> quadraticModificationList;

    size_t transferType; /*!Type of transfer function to use in defining covectors for modifications. 0 = dark matter, 1 = baryons*/

    ModificationManager(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext_,
                        const cosmology::CosmologicalParameters<T> &cosmology_,
                        std::vector<std::shared_ptr<fields::OutputField<DataType>>> &outputFields_) :
        outputFields(outputFields_), underlying(multiLevelContext_), cosmology(cosmology_) {
            this->transferType = 0; // Always modifying dark matter, never baryons.
        }

    //! Calculate existing value of the quantity defined by name
    template<typename ... Args>
    T calculateCurrentValueByName(std::string name_, Args &&... args) {

      std::shared_ptr<Modification<DataType, T>> modification = getModificationFromName(name_,
                                                                                        std::forward<Args>(args)...);
      outputFields[transferType]->toFourier();
      T value = modification->calculateCurrentValue( *(outputFields[transferType]) );
      return value;
    }

    /*!
        \param type_ Modification can be relative to existing value or absolute
        \param target_ Absolute target or factor by which the existing will be multiplied
      */
    template<typename ... Args>
    void addModificationToList(std::string name_, std::string type_, T target_, Args &&... args) {

      std::shared_ptr<Modification<DataType, T>> modification = getModificationFromName(name_,
                                                                                        std::forward<Args>(args)...);

      bool relative = isRelative(type_);
      T target = target_;

      if (relative) {
        T value = modification->calculateCurrentValue( *(outputFields[transferType]) );
        target *= value;
      }

      modification->setTarget(target);

      if (modification->getOrder() == 1) {
        linearModificationList.push_back(std::dynamic_pointer_cast<LinearModification<DataType, T>>(modification));
      } else if (modification->getOrder() == 2) {
        quadraticModificationList.push_back(
            std::dynamic_pointer_cast<QuadraticModification<DataType, T>>(modification));
      } else {
        throw std::runtime_error(" Could not add modification to list");
      }
    }

    bool hasModifications(){
      return !this->linearModificationList.empty() || !this->quadraticModificationList.empty();
    }

    //! Construct the modified field with all modifications present in the modification list
    void applyModifications() {

      std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas;
      std::vector<T> linear_targets;
      T pre_modif_chi2_from_field ;
      T post_modif_chi2_from_field ;

      pre_modif_chi2_from_field = outputFields[transferType]->getChi2();
      std::cerr << "BEFORE modifications chi^2 = " << pre_modif_chi2_from_field << std::endl;


      // Extract A, b from modification list
      for (size_t i = 0; i < linearModificationList.size(); i++) {
        alphas.push_back(linearModificationList[i]->getCovector());
        linear_targets.push_back(linearModificationList[i]->getTarget());
        std::cerr << "linear_targets[" << i << "] = " << linear_targets[i] << std::endl;
      }

      // Apply all linear modifications
      orthonormaliseModifications(alphas, linear_targets);
      std::cerr << "Delta chi^2 from linear modifications = " << getDeltaChi2FromLinearModifs( *(outputFields[transferType]) , alphas, linear_targets) << std::endl;
      applyLinearModif( *(outputFields[transferType]) , alphas, linear_targets);


      // Apply the joint linear and quadratic in an iterative procedure
      // TODO DeltaChi2 from quadratic could be calculated from the iterations and added here.
      applyLinQuadModif(alphas);

      post_modif_chi2_from_field = outputFields[transferType]->getChi2();
      std::cerr << "AFTER  modifications chi^2 = " << post_modif_chi2_from_field << std::endl;
      std::cerr << "         Total delta chi^2 = " << post_modif_chi2_from_field - pre_modif_chi2_from_field << std::endl;

      // Propagate them modification to the other fields if present:
      for(size_t i = 1;i < this->outputFields.size();i++) {
        propagateModifications(*(outputFields[i]));
      }
    }


    //!\brief Propagate modifications of this field to another field:
    void propagateModifications(fields::OutputField<DataType>& otherField)
    {
        //Must be in Fourier space to do this:
        outputFields[transferType]->toFourier();
        otherField.toFourier();

        T pre_modif_chi2_from_field = otherField.getChi2();
        std::cerr << "Propagating modification to field " << otherField.transferType << std::endl
                  << ". BEFORE modifications chi^2 = " << pre_modif_chi2_from_field << std::endl;

        //First, copy across the modified field:
        otherField.copyData( *(outputFields[transferType]) );

        //Divide by modifiedField power spectrum and apply outputField's power spectrum:
        otherField.applyTransferRatio(outputFields[transferType]->transferType);
        std::cerr << "otherField.transferType = " << otherField.transferType << std::endl;

        T post_modif_chi2_from_field = otherField.getChi2();
        std::cerr << "AFTER  modifications chi^2 = " << post_modif_chi2_from_field << std::endl;
        std::cerr << "         Total delta chi^2 = " << post_modif_chi2_from_field - pre_modif_chi2_from_field << std::endl;
    }

    void clearModifications() {
      std::cout << "Clearing modification list" << std::endl;
      linearModificationList.clear();
      quadraticModificationList.clear();
    }


  private:
    template<typename ... Args>
    std::shared_ptr<Modification<DataType, T>> getModificationFromName(std::string name_, Args &&... args) {
      try {
        auto modification = getLinearModificationFromName(name_);
        return modification;
      } catch (UnknownModificationException &e) {
        // If modification is unknown, it might be quadratic so swallow exception for now.
        try {
          auto modification = getQuadraticModificationFromName(name_, std::forward<Args>(args)...);
          return modification;
        } catch (UnknownModificationException &e2) {
          throw e2;
        }
      }

    }


    // Note - these functions would have to pass transferType if we wanted to modify baryons too. Currently
    // the code doesn't support this, so it isn't an issue, but it would have to be taken into account if
    // baryon modifications were implemented.
    std::shared_ptr<LinearModification<DataType, T>> getLinearModificationFromName(std::string name_) {
      if ((strcasecmp(name_.c_str(), "overdensity") == 0)) {
        return make_shared<OverdensityModification<DataType, T>>(underlying, cosmology);
      } else if ((strcasecmp(name_.c_str(), "potential") == 0)) {
        return make_shared<PotentialModification<DataType, T>>(underlying, cosmology);
      } else if ((strcasecmp(name_.c_str(), "vx") == 0)) {
        return make_shared<VelocityModification<DataType, T>>(underlying, cosmology, 0);
      } else if ((strcasecmp(name_.c_str(), "vy") == 0)) {
        return make_shared<VelocityModification<DataType, T>>(underlying, cosmology, 1);
      } else if ((strcasecmp(name_.c_str(), "vz") == 0)) {
        return make_shared<VelocityModification<DataType, T>>(underlying, cosmology, 2);
      }else {
        throw UnknownModificationException(name_ + " " + "is an unknown modification name");
      }
    }

    template<typename ... Args>
    std::shared_ptr<QuadraticModification<DataType, T>>
    getQuadraticModificationFromName(std::string name_, Args &&... args) {
      if ((strcasecmp(name_.c_str(), "variance") == 0)) {
        return make_shared<FilteredVarianceModification<DataType, T>>(underlying, cosmology,
                                                                      std::forward<Args>(args) ...);
      } else {
        throw UnknownModificationException(name_ + " " + "is an unknown modification name");
      }
    }


    /*!
     * Linear modifications are applied by orthonormalisation and adding the
     * correction term. See Roth et al 2016 for details
     */
    void applyLinearModif(fields::OutputField<DataType> &field,
                          std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas,
                          std::vector<T> &targets) {

      std::vector<T> existing_values;
      for (size_t i = 0; i < linearModificationList.size(); i++) {
        existing_values.push_back(linearModificationList[i]->calculateCurrentValue(field));
      }

      for (size_t i = 0; i < alphas.size(); i++) {
        auto &alpha_i = *(alphas[i]);
        auto dval_i = targets[i] - existing_values[i];

        // Constraints are essentially covectors. Convert to a vector, with correct weighting/filtering.
        alpha_i.convertToVector();

        alpha_i.toFourier(); // almost certainly already is in Fourier space, but just to be safe
        field.addScaled(alpha_i, dval_i);
        alpha_i.convertToCovector();
      }
    }

    void applyLinQuadModif(std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas) {

      // If quadratic are independent, we can apply them one by one in a loop.
      assert(areQuadraticIndependent());

      size_t numberQuadraticModifs = quadraticModificationList.size();
      for (size_t i = 0; i < numberQuadraticModifs; i++) {
        auto modif_i = quadraticModificationList[i];
        int init_n_steps = modif_i->getInitNumberSteps();

        // Try the procedure on a test field and deduce the correct number of steps
        auto test_field = fields::OutputField<DataType>( *(outputFields[transferType]) );
        performIterations(test_field, alphas, modif_i, init_n_steps);
        int n_steps = calculateCorrectNumberSteps(test_field, modif_i, init_n_steps);

        // Perform procedure on real output
        if (n_steps > init_n_steps) {
          std::cout << n_steps << " steps are required for the quadratic algorithm " << std::endl;
          performIterations( *(outputFields[transferType]) , alphas, modif_i, n_steps);
        } else {
          std::cout << "No need to do more steps to achieve target precision" << std::endl;
          performIterations( *(outputFields[transferType]) , alphas, modif_i, init_n_steps);
        }


      }
    }

    //! Executes n_steps iterations of linear and quadratic modifications
    void performIterations(fields::OutputField<DataType> &field,
                           std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas,
                           std::shared_ptr<QuadraticModification<DataType, T>> quad_modif, int n_steps) {

      T overall_quad_target = quad_modif->getTarget();
      T starting_quad_value = quad_modif->calculateCurrentValue(field);

      std::vector<T> quad_targets = tools::linspace(starting_quad_value, overall_quad_target, n_steps);

      for (int i = 0; i < (n_steps); i++) {

        T current_value = quad_modif->calculateCurrentValue(field);
        auto pushedField = quad_modif->pushMultiLevelFieldThroughMatrix(field);
        pushedField->toFourier();

        T norm = sqrt(pushedField->innerProduct(*pushedField).real());
        addToOrthonormalFamily(alphas, pushedField);

        //Apply quad step
        T multiplier = 0.5 * (quad_targets[i + 1] - current_value) /
                       norm; //One sqrt factor inside the orthonormalise method and one more here.
        pushedField->convertToVector();
        field.addScaled(*pushedField, multiplier);
      }

    }

    int calculateCorrectNumberSteps(const fields::OutputField<DataType> &field,
                                    std::shared_ptr<QuadraticModification<DataType, T>> modif, int previous_n_steps) {

      T achieved_precision = std::abs(modif->calculateCurrentValue(field) - modif->getTarget());
      T target_precision = modif->getTarget() * modif->getTargetPrecision();

      int n_steps = previous_n_steps * (int) std::ceil(std::sqrt(achieved_precision / target_precision));
      return n_steps;
    }

    //! Graam-Schmidt procedure to orthonormalise the modification covectors
    void orthonormaliseModifications(std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas,
                                     std::vector<T> &targets) {

      using namespace tools::numerics;
      size_t n = alphas.size();

      // Calculate the inner products
      for (size_t i = 0; i < n; i++) {
        auto &alpha_i = *(alphas[i]);
        for (size_t j = 0; j < i; j++) {
          auto &alpha_j = *(alphas[j]);
          T result = alpha_i.innerProduct(alpha_j).real();

          alpha_i.addScaled(alpha_j, -result);

          // update constraining value
          targets[i] -= result * targets[j];
        }

        // normalize
        T norm = sqrt(alpha_i.innerProduct(alpha_i).real());
        alpha_i /= norm;
        targets[i] /= norm;
      }
    }

    //! Orthonormalise a covector with respect to an already orthonormal family
    void addToOrthonormalFamily(std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas,
                                std::shared_ptr<fields::ConstraintField<DataType>> alpha) {

      using namespace tools::numerics;
      size_t n = alphas.size();

      // Calculate the inner products between the new vector and the existing family
      for (size_t i = 0; i < n; i++) {
        auto &alpha_i = *(alphas[i]);
        T result = alpha->innerProduct(alpha_i).real();
        alpha->addScaled(alpha_i, -result);
      }

      // normalize
      T norm = sqrt(alpha->innerProduct(*alpha).real());
      (*alpha) /= norm;
    }

    bool isRelative(std::string type) {
      bool relative = false;
      if (strcasecmp(type.c_str(), "relative") == 0) {
        relative = true;
      } else if (strcasecmp(type.c_str(), "absolute") != 0) {
        throw std::runtime_error("Modification type must be either 'relative' or 'absolute'");
      }
      return relative;
    }

    //! Calculates delta chi2 from orthonormalised linear modifications (Eq 12 Roth et al 2016)
    T getDeltaChi2FromLinearModifs(fields::OutputField<DataType> &field,
                              std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> alphas,
                              std::vector<T> &targets){

      T target_norm = T(0.0);
      T initial_values_norm = T(0.0);

      for (size_t i = 0; i < alphas.size(); i++) {
        auto &alpha_i = *(alphas[i]);
        initial_values_norm += std::pow(std::abs(alpha_i.innerProduct(field).real()), 2);
        target_norm += std::pow(std::abs(targets[i]), 2);
      }

      return target_norm - initial_values_norm;
    }

    //! Check quadratic modifications are independent by calculating delta * Q1 * C0 * Q2 * delta
    bool areQuadraticIndependent(){
      //TODO Condition here is necessary but not sufficient to ensure independence.
      // Not sure if we ought to keep this anyway and leave independence as a user problem

      size_t n = quadraticModificationList.size();
      bool indep = true;

      for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < i; j++){
          auto pushed_i = quadraticModificationList[i]->pushMultiLevelFieldThroughMatrix( *(outputFields[transferType]) );
          auto pushed_j = quadraticModificationList[j]->pushMultiLevelFieldThroughMatrix( *(outputFields[transferType]) );

          auto ortho = pushed_i->innerProduct(*pushed_j).real();

          if(ortho > 0.001){
            indep = false ;
          }
        }
      }
      return indep;
    }
  };
}
#endif
