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
  class ZeldovichParticleGenerator : public ParticleGenerator<GridDataType> {
  protected:
    using TField = fields::Field<GridDataType, T>;
    using TRealField = fields::Field<T,T>;

    const cosmology::CosmologicalParameters<T> &cosmology;
    TField &linearOverdensityField;
    using ParticleGenerator<GridDataType>::grid;
    T velocityToOffsetRatio, boxMass;

    // The grid offsets after Zeldovich approximation is applied
    // (nullptr before that):
    std::shared_ptr<TRealField> pOff_x;
    std::shared_ptr<TRealField> pOff_y;
    std::shared_ptr<TRealField> pOff_z;



    void calculateVelocityToOffsetRatio() {

      velocityToOffsetRatio = 1. * 100. * sqrt(
        cosmology.OmegaM0 / cosmology.scalefactor / cosmology.scalefactor / cosmology.scalefactor +
        cosmology.OmegaLambda0) * sqrt(
        cosmology.scalefactor);

      //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
      //TODO: hardcoded value of f=1 is inaccurate - should be a function of omega

    }

    void calculateSimulationMass() {
      // The mass of the entire simulation:
      // sqrt(3*(100 h km/s/Mpc)^2/(8 pi G))
      //
      // Gadget units are used internally, so express as Msol h
      boxMass = 27.744948 * cosmology.OmegaM0 * powf(grid.simsize, 3.0);
    }

    virtual void calculateOffsetFields() {
      // TODO: refactorise this horrible long method

      size_t size = grid.size;
      tools::progress::ProgressBar pb("zeldovich", size * 2);

      // get a reference to the density field in fourier space
      linearOverdensityField.toFourier();

      // copy three times to start assembling the vx, vy, vz fields
      auto offsetX = std::make_shared<TField>(const_cast<grids::Grid<T> &>(grid));
      auto offsetY = std::make_shared<TField>(const_cast<grids::Grid<T> &>(grid));
      auto offsetZ = std::make_shared<TField>(const_cast<grids::Grid<T> &>(grid));

      const T kw = 2. * M_PI / grid.boxsize;
      const int nyquist = tools::numerics::fourier::getNyquistModeThatMustBeReal(grid);

      tools::numerics::fourier::applyTransformationInFourierBasis<T>(linearOverdensityField,
      [kw, nyquist](complex<T> inputVal, int iix, int iiy, int iiz) -> std::tuple<complex<T>, complex<T>, complex<T>> {
        complex<T> result_x;
        T kfft = (iix * iix + iiy * iiy + iiz * iiz);
        result_x.real(-inputVal.imag()/(kfft*kw));
        result_x.imag(inputVal.real()/(kfft*kw));
        complex<T> result_y(result_x);
        complex<T> result_z(result_x);

        result_x*=iix;
        result_y*=iiy;
        result_z*=iiz;

        // derivative at nyquist frequency is not defined; set it to zero
        // potential is also undefined at (0,0,0); set that mode to zero too
        if(abs(iix)==nyquist || kfft==0)
          result_x=0;
        if(abs(iiy)==nyquist || kfft==0)
          result_y=0;
        if(abs(iiz)==nyquist || kfft==0)
          result_z=0;

        return std::make_tuple(result_x, result_y, result_z);
      },
      *offsetX, *offsetY, *offsetZ);

      offsetX->toReal();
      offsetY->toReal();
      offsetZ->toReal();

      // if the implementation is complex, at this point get real copies of the fields
      // otherwise just use the offset field we have already constructed
      this->pOff_x = getRealPart(*offsetX);
      this->pOff_y = getRealPart(*offsetY);
      this->pOff_z = getRealPart(*offsetZ);

    }

  public:


    ZeldovichParticleGenerator(TField &linearOverdensityField,
                               const cosmology::CosmologicalParameters<T> &cosmology) :

    ParticleGenerator<GridDataType>(linearOverdensityField.getGrid()),
    cosmology(cosmology),
    linearOverdensityField(linearOverdensityField)
    {
      recalculate();
    }

    void recalculate() override {
      calculateVelocityToOffsetRatio();
      calculateSimulationMass();
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

    virtual T getMass(const grids::Grid<T> & onGrid) const override {
      return boxMass*onGrid.cellMassFrac;
    }

    virtual T getEps(const grids::Grid<T> & onGrid) const override  {
      return onGrid.dx * onGrid.cellSofteningScale * 0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    auto getOffsetFields() {
      return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    auto getOffsetFields() const {
      return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    virtual particle::Particle<T> getParticleNoOffset(const grids::Grid<T> &onGrid, size_t id) const override {
      particle::Particle<T> particle;

      assert(!pOff_x->isFourier());
      assert(!pOff_y->isFourier());
      assert(!pOff_z->isFourier());

      particle.pos.x = onGrid.getFieldAt(id, *pOff_x);
      particle.pos.y = onGrid.getFieldAt(id, *pOff_y);
      particle.pos.z = onGrid.getFieldAt(id, *pOff_z);

      particle.vel = particle.pos*velocityToOffsetRatio;

      particle.mass = getMass(onGrid);
      particle.soft = getEps(onGrid);

      return particle;
    }

    virtual particle::Particle<T> getParticleNoWrap(const grids::Grid<T> &onGrid, size_t id) const override {
      auto particle = getParticleNoOffset(onGrid, id);
      auto centroid = onGrid.getCellCentroid(id);
      particle.pos+=centroid;
      return particle;
    }



  };
}

#endif //IC_ZELDOVICH_HPP
