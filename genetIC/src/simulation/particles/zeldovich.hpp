//
// Created by Andrew Pontzen on 20/12/2016.
//

#ifndef IC_ZELDOVICH_HPP
#define IC_ZELDOVICH_HPP

#include <complex>
#include <src/io/numpy.hpp>
#include "src/tools/progress/progress.hpp"
#include "src/simulation/field/field.hpp"
#include "src/simulation/particles/generator.hpp"

namespace cosmology {
  template<typename T>
  struct CosmologicalParameters;
}

namespace particle {

  template<typename T>
  class Particle;

  template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
  class ZeldovichParticleGenerator;

  /*! \class ZeldovichParticleEvaluator
      \brief Class to evaluate particles on a grid using the Zeldovich method
  */
  template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
  class ZeldovichParticleEvaluator : public ParticleEvaluator<GridDataType> {
  protected:
    using EvaluatorType = std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>;
    using GridType = grids::Grid<T>;
    using TField = fields::Field<GridDataType, T>;
    EvaluatorType pOffsetXEvaluator; //!< Evaluator for the x offsets
    EvaluatorType pOffsetYEvaluator; //!< Evaluator for the y offsets
    EvaluatorType pOffsetZEvaluator; //!< Evaluator for the z offsets

    const cosmology::CosmologicalParameters<T> &cosmology; //!< Cosmological parameters

    T velocityToOffsetRatio; //!< velocity offset ratio, for converting position offsets into velocity offsets
    T boxMass; //!< Total matter mass in simulation
    T epsNorm; //!< Normalisation for the cell softening scale for generated particles

    std::shared_ptr<const GridType> onGrid; //!< Pointer to the grid for this evaluator

    //! Computes the ratio between position and velocity offsets for the given cosmological parameters
    void calculateVelocityToOffsetRatio() {
      velocityToOffsetRatio = cosmology::zeldovichVelocityToOffsetRatio(cosmology);
    }

    //! Computes the total mass in the simulation
    void calculateSimulationMass() {
      // The mass of the entire simulation:
      // sqrt(3*(100 h km/s/Mpc)^2/(8 pi G))
      //
      // Gadget units are used internally, so express as Msol h
      // To be precise, this pre-factor is the critical density in units of (10^10*Msol//)/(Mpc/h)^3, where
      // Msol is the mass of the sun ('solar mass'). Needed because volumes are expected in units of Mpc/h
      boxMass = 27.744948 * cosmology.OmegaM0 * pow(onGrid->periodicDomainSize, 3.0);
    }

  public:
    /*! \brief Constructor that creates the evaluator based on the evaluators for each of the x,y,z positions, and the specified grid and cell softening scale
        \param evalOffX - evaluator for the x offsets
        \param evalOffX - evaluator for the y offsets
        \param evalOffX - evaluator for the z offsets
        \param grid - underlying grid
        \param cosmology - cosmological parameters
        \param epsNorm_ - prefactor used to define cell softening scale
    */
    ZeldovichParticleEvaluator(EvaluatorType evalOffX, EvaluatorType evalOffY, EvaluatorType evalOffZ,
                               const GridType &grid, const cosmology::CosmologicalParameters<T> &cosmology,
                               T epsNorm_ = 0.01075) // Using default value (arbitrary to coincide with normal UW resolution)
      : ParticleEvaluator<GridDataType>(grid), cosmology(cosmology) {
      pOffsetXEvaluator = evalOffX;
      pOffsetYEvaluator = evalOffY;
      pOffsetZEvaluator = evalOffZ;
      onGrid = grid.shared_from_this();
      calculateSimulationMass();
      calculateVelocityToOffsetRatio();
      epsNorm = epsNorm_;
    }

    //! Evaluates the particle at cell id, using the x,y,z offset evaluators, and the appropriate velocity offset ratio
    virtual particle::Particle<T> getParticleNoOffset(size_t id) const override {

      particle::Particle<T> particle;

      particle.pos.x = tools::datatypes::real_part_if_complex((*pOffsetXEvaluator)[id]);
      particle.pos.y = tools::datatypes::real_part_if_complex((*pOffsetYEvaluator)[id]);
      particle.pos.z = tools::datatypes::real_part_if_complex((*pOffsetZEvaluator)[id]);

      particle.vel = particle.pos * velocityToOffsetRatio;

      particle.mass = getMass();
      particle.soft = getEps();

      return particle;
    }

    //! Evaluates particle without wrapping
    virtual particle::Particle<T> getParticleNoWrap(size_t id) const override {
      auto particle = getParticleNoOffset(id);
      auto centroid = onGrid->getCentroidFromIndex(id);
      particle.pos += centroid;
      return particle;
    }

    //! Gets the mass for a single particle
    virtual T getMass() const override {
      return boxMass * onGrid->cellMassFrac;
    }

    //! Gets the cell softening scale for a single particle
    virtual T getEps() const override {
      return onGrid->cellSize * onGrid->cellSofteningScale * epsNorm;
    }

  };


  /*! \class ZeldovichParticleGenerator
      \brief Class to generate particles using the Zeldovich approximation
  */
  template<typename GridDataType, typename T>
  class ZeldovichParticleGenerator : public ParticleGenerator<GridDataType> {
  protected:
    using TField = fields::Field<GridDataType, T>;
    using TRealField = fields::Field<T, T>;

    friend class ZeldovichParticleEvaluator<GridDataType, T>;

    TField &linearOverdensityField; //!< Overdensity field used to generate particles on this grid
    using ParticleGenerator<GridDataType>::grid;


    // The grid offsets after Zeldovich approximation is applied
    // (nullptr before that):
    std::shared_ptr<TField> pOff_x; //!< Offset field for x positions
    std::shared_ptr<TField> pOff_y; //!< Offset field for y positions
    std::shared_ptr<TField> pOff_z; //!< Offset field for z positions


    //! Calculates the offset fields, which will be used to actually compute the position and velocity offsets by the Zeldovich evaluator
    virtual void calculateOffsetFields() {
      size_t size = grid.size;

      const T nyquist = tools::numerics::fourier::getNyquistModeThatMustBeReal(grid) * grid.getFourierKmin();

      auto zeldovichOffsetFields = linearOverdensityField.generateNewFourierFields(
        [nyquist](complex<T> inputVal, T kx, T ky, T kz) -> std::tuple<complex<T>, complex<T>, complex<T>> {
          complex<T> result_x;
          T kfft = kx * kx + ky * ky + kz * kz; // k^2

          // Computes i*inputVal/k^2:
          result_x.real(-inputVal.imag() / (kfft));
          result_x.imag(inputVal.real() / (kfft));
          complex<T> result_y(result_x);
          complex<T> result_z(result_x);

          result_x *= kx; // i*kx*inputVal/k^2
          result_y *= ky;
          result_z *= kz;

          // derivative at nyquist frequency is not defined; set it to zero
          // potential is also undefined at (0,0,0); set that mode to zero too
          if (kx == nyquist || kfft == 0)
            result_x = 0;
          if (ky == nyquist || kfft == 0)
            result_y = 0;
          if (kz == nyquist || kfft == 0)
            result_z = 0;

          return std::make_tuple(result_x, result_y, result_z);
        });

      std::tie(this->pOff_x, this->pOff_y, this->pOff_z) = zeldovichOffsetFields;
      this->pOff_x->toReal();
      this->pOff_y->toReal();
      this->pOff_z->toReal();


    }

  public:


    //! Constructor from a given overdensity field
    ZeldovichParticleGenerator(TField &linearOverdensityField) :
      ParticleGenerator<GridDataType>(linearOverdensityField.getGrid()),
      linearOverdensityField(linearOverdensityField) {
      recalculate();
    }

    //! Computes the offset fields using the Zeldovich approximation
    void recalculate() override {
      calculateOffsetFields();
    }

    //! Adds offset fields from another generateo defined on a different grid to the ones generated by this generator
    void addFieldFromDifferentGrid(const ZeldovichParticleGenerator &source) {
      pOff_x->addFieldFromDifferentGrid(*source.pOff_x);
      pOff_y->addFieldFromDifferentGrid(*source.pOff_y);
      pOff_z->addFieldFromDifferentGrid(*source.pOff_z);
    }

    //! Adds offset fields from a different grid, but applies a filter to them first
    void addFieldFromDifferentGridWithFilter(ZeldovichParticleGenerator &source, const filters::Filter<T> &filter) {
      pOff_x->addFieldFromDifferentGridWithFilter(*source.pOff_x, filter);
      pOff_y->addFieldFromDifferentGridWithFilter(*source.pOff_y, filter);
      pOff_z->addFieldFromDifferentGridWithFilter(*source.pOff_z, filter);
    }

    //! Applies a filter to the offset fields
    void applyFilter(const filters::Filter<T> &filter) {
      pOff_x->applyFilter(filter);
      pOff_y->applyFilter(filter);
      pOff_z->applyFilter(filter);
    }

    //! Fourier transforms the offset fields to real space (they are initially defined in Fourier space)
    void toReal() {
      pOff_x->toReal();
      pOff_y->toReal();
      pOff_z->toReal();
    }

    //! Returns a tuple with the created offset fields
    std::vector<std::shared_ptr<fields::Field<GridDataType>>> getGeneratedFields() override {
      return {pOff_x, pOff_y, pOff_z};
    }


  };
}

#endif //IC_ZELDOVICH_HPP
