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

  template<typename GridDataType, typename T=tools::datatypes::strip_complex<GridDataType>>
  class ZeldovichParticleEvaluator : public ParticleEvaluator<GridDataType> {
  protected:
    using EvaluatorType = std::shared_ptr<fields::EvaluatorBase<GridDataType, T>>;
    using GridType = grids::Grid<T>;
    using TField = fields::Field<GridDataType, T>;
    EvaluatorType pOffsetXEvaluator;
    EvaluatorType pOffsetYEvaluator;
    EvaluatorType pOffsetZEvaluator;

    const cosmology::CosmologicalParameters<T> &cosmology;

    T velocityToOffsetRatio, boxMass, epsNorm;

    std::shared_ptr<const GridType> onGrid;

    void calculateVelocityToOffsetRatio() {

    //Just to make this easier to read (should be optimised away by the compiler anyway):
      T f = 1.0;
      T a = cosmology.scalefactor;
      T Om = cosmology.OmegaM0;
      T Ol = cosmology.OmegaLambda0;
      //According to Carrol and Press (1992) should use:
      //f = powf(( Om/(a*a*a) )/( Om/(a*a*a) + (1.0 - Om - Ol)/(a*a) + Ol ),4.0/7.0);

      velocityToOffsetRatio = f * 100. * sqrt( ( Om / (a*a*a) ) +
          // + (1.0 - cosmology.OmegaM0 - cosmology.OmegaLambda0)/(a*a) //We use the curvature term elsewhere, so why not here?
          Ol) * sqrt(a);

      //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
      //TODO: hardcoded value of f=1 is inaccurate - should be a function of omega

    }

    void calculateSimulationMass() {
      // The mass of the entire simulation:
      // sqrt(3*(100 h km/s/Mpc)^2/(8 pi G))
      //
      // Gadget units are used internally, so express as Msol h
      //To be precise, this pre-factor is the critical density in units of (10^10*Msol//)/(Mpc/h)^3, where
      //Msol is the mass of the sun ('solar mass'). Needed because volumes are expected in units of Mpc/h
      //But - shouldn't this total mass actually include the mass of baryons?
      boxMass = 27.744948 * cosmology.OmegaM0 * powf(onGrid->periodicDomainSize, 3.0);
    }

  public:
    ZeldovichParticleEvaluator(EvaluatorType evalOffX, EvaluatorType evalOffY, EvaluatorType evalOffZ,
                               const GridType &grid, const cosmology::CosmologicalParameters<T> &cosmology,
                               T epsNorm_ = 0.01075)//Using default value (arbitrary to coincide with normal UW resolution)
        : ParticleEvaluator<GridDataType>(grid), cosmology(cosmology) {
      pOffsetXEvaluator = evalOffX;
      pOffsetYEvaluator = evalOffY;
      pOffsetZEvaluator = evalOffZ;
      onGrid = grid.shared_from_this();
      calculateSimulationMass();
      calculateVelocityToOffsetRatio();
      epsNorm = epsNorm_;
    }

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

    virtual particle::Particle<T> getParticleNoWrap(size_t id) const override {
      auto particle = getParticleNoOffset(id);
      auto centroid = onGrid->getCentroidFromIndex(id);
      particle.pos += centroid;
      return particle;
    }

    virtual T getMass() const override {
      return boxMass * onGrid->cellMassFrac;
    }

    virtual T getEps() const override {
      return onGrid->cellSize * onGrid->cellSofteningScale * epsNorm;
             //0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

  };


  template<typename GridDataType, typename T>
  class ZeldovichParticleGenerator : public ParticleGenerator<GridDataType> {
  protected:
    using TField = fields::Field<GridDataType, T>;
    using TRealField = fields::Field<T, T>;

    friend class ZeldovichParticleEvaluator<GridDataType, T>;

    TField &linearOverdensityField;
    using ParticleGenerator<GridDataType>::grid;


    // The grid offsets after Zeldovich approximation is applied
    // (nullptr before that):
    std::shared_ptr<TField> pOff_x;
    std::shared_ptr<TField> pOff_y;
    std::shared_ptr<TField> pOff_z;


    virtual void calculateOffsetFields() {
      size_t size = grid.size;
      tools::progress::ProgressBar pb("zeldovich", size * 2);

      const T nyquist = tools::numerics::fourier::getNyquistModeThatMustBeReal(grid) * grid.getFourierKmin();

      auto zeldovichOffsetFields = linearOverdensityField.generateNewFourierFields(
          [nyquist](complex<T> inputVal, T kx, T ky, T kz) -> std::tuple<complex<T>, complex<T>, complex<T>> {
            complex<T> result_x;
            T kfft = kx * kx + ky * ky + kz * kz;

            result_x.real(-inputVal.imag() / (kfft));
            result_x.imag(inputVal.real() / (kfft));
            complex<T> result_y(result_x);
            complex<T> result_z(result_x);

            result_x *= kx;
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


    ZeldovichParticleGenerator(TField &linearOverdensityField) :
        ParticleGenerator<GridDataType>(linearOverdensityField.getGrid()),
        linearOverdensityField(linearOverdensityField) {
      recalculate();
    }

    void recalculate() override {
      calculateOffsetFields();
    }

    void addFieldFromDifferentGrid(const ZeldovichParticleGenerator &source) {
      pOff_x->addFieldFromDifferentGrid(*source.pOff_x);
      pOff_y->addFieldFromDifferentGrid(*source.pOff_y);
      pOff_z->addFieldFromDifferentGrid(*source.pOff_z);
    }

    void addFieldFromDifferentGridWithFilter(ZeldovichParticleGenerator &source, const filters::Filter<T> &filter) {
      pOff_x->addFieldFromDifferentGridWithFilter(*source.pOff_x, filter);
      pOff_y->addFieldFromDifferentGridWithFilter(*source.pOff_y, filter);
      pOff_z->addFieldFromDifferentGridWithFilter(*source.pOff_z, filter);
    }

    void applyFilter(const filters::Filter<T> &filter) {
      pOff_x->applyFilter(filter);
      pOff_y->applyFilter(filter);
      pOff_z->applyFilter(filter);
    }

    void toReal() {
      pOff_x->toReal();
      pOff_y->toReal();
      pOff_z->toReal();
    }

    std::vector<std::shared_ptr<fields::Field<GridDataType>>> getGeneratedFields() override {
      return {pOff_x, pOff_y, pOff_z};
    }


  };
}

#endif //IC_ZELDOVICH_HPP
