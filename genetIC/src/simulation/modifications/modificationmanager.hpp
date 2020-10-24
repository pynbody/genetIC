#ifndef IC_MODIFICATIONMANAGER_HPP
#define IC_MODIFICATIONMANAGER_HPP

#include "src/simulation/modifications/linearmodification.hpp"
#include "src/simulation/modifications/quadraticmodification.hpp"
#include "src/tools/logging.hpp"
#include <string>

//! Deals with the creation of genetically modified fields
namespace modifications {

  //! Keep track of modifications to be applied. Construct the modified field from these modifications.
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  class ModificationManager {

  protected:
    std::shared_ptr<fields::OutputField<DataType>> outputField;      //!< The field on which modifications are being made
    const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext;          //!< Grid context in which modifications take place
    const cosmology::CosmologicalParameters<T> &cosmology;                          //!< Cosmology context in which modifications take place
    std::vector<std::shared_ptr<LinearModification<DataType, T>>> linearModificationList;  //!< Modifications to be applied
    std::vector<std::shared_ptr<QuadraticModification<DataType, T>>> quadraticModificationList; //!< List of quadratic modifications to be applied

  public:
    //! \brief Constructor which accepts a multi-level context, cosmological parameters, and the fields to modify
    /*! \param multiLevelContext_ - reference to the multi-level context object
        \param cosmology_ - stores cosmological parameters
        \param outputFields_ - vector of shared pointers to the output fields. Only the dark matter is modified, but we need to propagate the modifications to the others.
    */
    ModificationManager(multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext_,
                        const cosmology::CosmologicalParameters<T> &cosmology_,
                        const std::shared_ptr<fields::OutputField<DataType>> &outputField_) :
      outputField(outputField_), multiLevelContext(multiLevelContext_), cosmology(cosmology_) {

    }

    //! Calculate existing value of the quantity defined by name
    template<typename ... Args>
    T calculateCurrentValueByName(std::string name_, Args &&... args) {

      std::shared_ptr<Modification<DataType, T>> modification = getModificationFromName(name_,
                                                                                        std::forward<Args>(args)...);
      outputField->toFourier();
      T value = modification->calculateCurrentValue(*outputField);
      return value;
    }

    //! Specify the field to which modifications are to be made
    void bindToField(const std::shared_ptr<fields::OutputField<DataType>> &field) {
      this->outputField = field;
    }

    /*! \brief Adds the specified modification to the list of modifications to be applied, if it exists.

        \param name_ - String naming the modification required.
        \param type_ - Modification can be relative to existing value or absolute.
        \param target_ - Absolute target or factor by which the existing will be multiplied.
      */
    template<typename ... Args>
    void addModificationToList(std::string name_, std::string type_, T target_, Args &&... args) {

      std::shared_ptr<Modification<DataType, T>> modification = getModificationFromName(name_,
                                                                                        std::forward<Args>(args)...);

      bool relative = isRelative(type_);
      T target = target_;

      if (relative) {
        T value = modification->calculateCurrentValue(*outputField);
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

    //! Returns true if modifications have been added to the list
    bool hasModifications() {
      return !this->linearModificationList.empty() || !this->quadraticModificationList.empty();
    }

    //! Construct the modified field with all modifications present in the modification list, then propagate it to the other fields.
    void applyModifications() {

      std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> modificationCovectors;
      std::vector<T> linearTargetValues;
      T pre_modif_chi2_from_field;
      T post_modif_chi2_from_field;

      pre_modif_chi2_from_field = outputField->getChi2();


      // Extract A, b from modification list
      for (size_t i = 0; i < linearModificationList.size(); i++) {
        modificationCovectors.push_back(linearModificationList[i]->getCovector(outputField->getTransferType()));
        linearTargetValues.push_back(linearModificationList[i]->getTarget());
      }

      // Apply all linear modifications
      logging::entry() << std::endl << "Applying modifications" << std::endl;
      orthonormaliseModifications(modificationCovectors, linearTargetValues);
#ifdef DEBUG_INFO
      logging::entry() << "ESTIMATED delta chi^2 from all linear modifications = "
                << getDeltaChi2FromLinearModifs(*outputField, modificationCovectors, linearTargetValues)
                << std::endl;
#endif

      applyLinearModif( modificationCovectors, linearTargetValues);
      applyLinQuadModif(modificationCovectors);

      post_modif_chi2_from_field = outputField->getChi2();
      logging::entry() << "   Post-modification chi^2 = " << post_modif_chi2_from_field << std::endl;
      size_t dof = this->multiLevelContext.getNumDof();
      logging::entry() << "  Modification Delta chi^2 = " << post_modif_chi2_from_field - pre_modif_chi2_from_field
                << std::endl;
      logging::entry() << "           d.o.f. in field = " << dof << std::endl;
      logging::entry() << std::endl;
    }


    //! Clear all modifications from the list of modifications to be applied.
    void clearModifications() {
      logging::entry() << "Clearing modification list" << std::endl;
      linearModificationList.clear();
      quadraticModificationList.clear();
    }


  private:
    //! Returns the modification a supplied string
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
    //! Returns the appropriate linear modification from the supplied string, if it exists
    std::shared_ptr<LinearModification<DataType, T>> getLinearModificationFromName(std::string name_) {
      if ((strcasecmp(name_.c_str(), "overdensity") == 0)) {
        return make_shared<OverdensityModification<DataType, T>>(multiLevelContext, cosmology);
      } else if ((strcasecmp(name_.c_str(), "potential") == 0)) {
        return make_shared<PotentialModification<DataType, T>>(multiLevelContext, cosmology);
      } else if ((strcasecmp(name_.c_str(), "vx") == 0)) {
        return make_shared<VelocityModification<DataType, T>>(multiLevelContext, cosmology, 0);
      } else if ((strcasecmp(name_.c_str(), "vy") == 0)) {
        return make_shared<VelocityModification<DataType, T>>(multiLevelContext, cosmology, 1);
      } else if ((strcasecmp(name_.c_str(), "vz") == 0)) {
        return make_shared<VelocityModification<DataType, T>>(multiLevelContext, cosmology, 2);
      } else {
        throw UnknownModificationException(name_ + " " + "is an unknown modification name");
      }
    }

    //! Returns the appropriate quadratic modification from the supplied string, if it exists
    template<typename ... Args>
    std::shared_ptr<QuadraticModification<DataType, T>>
    getQuadraticModificationFromName(std::string name_, Args &&... args) {
      if ((strcasecmp(name_.c_str(), "variance") == 0)) {
        return make_shared<FilteredVarianceModification<DataType, T>>(multiLevelContext, cosmology,
                                                                      std::forward<Args>(args) ...);
      } else {
        throw UnknownModificationException(name_ + " " + "is an unknown modification name");
      }
    }


    /*!
     * Linear modifications are applied by orthonormalisation and adding the
     * correction term. See Roth et al 2016 for details
     */
    void applyLinearModif(std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> orthonormalisedCovectors,
                          std::vector<T> &orthonormalisedTargetValues) {

      std::vector<T> existingValues;
      for (size_t i = 0; i < orthonormalisedCovectors.size(); i++) {
        existingValues.push_back(orthonormalisedCovectors[i]->innerProduct(*outputField).real());
      }

      for (size_t i = 0; i < orthonormalisedCovectors.size(); i++) {
        auto &alpha_i = *(orthonormalisedCovectors[i]);
        auto dval_i = orthonormalisedTargetValues[i] - existingValues[i];

        // Convert to a vector, with correct weighting/filtering, in preparation for adding to the output field
        alpha_i.convertToVector();

        alpha_i.toFourier(); // almost certainly already is in Fourier space, but just to be safe
        outputField->addScaled(alpha_i, dval_i);
      }
    }

    //! Apply quadratic modifications, while keeping values of the specified covectors fixed
    void applyLinQuadModif(std::vector<std::shared_ptr<fields::ConstraintField<DataType>>> orthonormalisedCovectors) {

      // If quadratic are independent, we can apply them one by one in a loop.
      assert(areQuadraticIndependent());

      size_t numberQuadraticModifs = quadraticModificationList.size();
      for (size_t i = 0; i < numberQuadraticModifs; i++) {
        auto modif_i = quadraticModificationList[i];
        int init_n_steps = modif_i->getInitNumberSteps();

        // Try the procedure on a test field and deduce the correct number of steps
        auto test_field = fields::OutputField<DataType>(*outputField);
        performIterations(test_field, orthonormalisedCovectors, modif_i, init_n_steps);
        int n_steps = calculateCorrectNumberSteps(test_field, modif_i, init_n_steps);

        // Perform procedure on real output
        if (n_steps > init_n_steps) {
          logging::entry() << n_steps << " steps are required for the quadratic algorithm " << std::endl;
          performIterations(*outputField, orthonormalisedCovectors, modif_i, n_steps);
        } else {
          logging::entry() << "No need to do more steps to achieve target precision" << std::endl;
          performIterations(*outputField, orthonormalisedCovectors, modif_i, init_n_steps);
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

    //! Compute number of steps needed to apply a quadratic modification.
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

    //! Returns true if the supplied string is "relative", false if "absolute", and an error if anything else.
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
                                   std::vector<T> &targets) {

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
    bool areQuadraticIndependent() {
      //TODO Condition here is necessary but not sufficient to ensure independence.
      // Not sure if we ought to keep this anyway and leave independence as a user problem

      size_t n = quadraticModificationList.size();
      bool indep = true;

      for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
          auto pushed_i = quadraticModificationList[i]->pushMultiLevelFieldThroughMatrix(*outputField);
          auto pushed_j = quadraticModificationList[j]->pushMultiLevelFieldThroughMatrix(*outputField);

          auto ortho = pushed_i->innerProduct(*pushed_j).real();

          if (ortho > 0.001) {
            indep = false;
          }
        }
      }
      return indep;
    }
  };
}
#endif
