//
// Created by Andrew Pontzen on 20/12/2016.
//

#ifndef IC_ZELDOVICH_HPP
#define IC_ZELDOVICH_HPP

#include <complex>
#include <src/numpy.hpp>
#include "../cosmo.hpp"
#include "../progress/progress.hpp"
#include "../field.hpp"
#include "generator.hpp"
#include "particle.hpp"

template<typename T>
struct CosmologicalParameters;

namespace particle {

  template<typename T>
  class Particle;

  template<typename GridDataType, typename T=strip_complex<GridDataType>>
  class ZeldovichParticleGenerator : public ParticleGenerator<GridDataType> {
  protected:
    using TField = Field<GridDataType, T>;
    using TRealField = Field<T,T>;

    T velocityToOffsetRatio, boxMass;
    TField &linearOverdensityField;
    const CosmologicalParameters<T> &cosmology;
    using ParticleGenerator<GridDataType>::grid;

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
      boxMass = 27.78 * cosmology.OmegaM0 * powf(grid.simsize, 3.0);
    }

    virtual void calculateOffsetFields() {
      // TODO: refactorise this horrible long method

      size_t size = grid.size;
      size_t size3 = grid.size3;
      progress::ProgressBar pb("zeldovich", size * 2);

      // get a reference to the density field in fourier space
      linearOverdensityField.toFourier();
      auto &pField_k = linearOverdensityField.getDataVector();

      /*
      RE-FFT:
      TEST:(-0.470325,-0.0842186)
(-0.417907,-0.314688)
(-0.299121,0.507641)
(-0.223375,0.369053)
(0.246137,-0.391641)

      ORIG-FFT:
      TEST:(-0.470325,-0.0842186)
(-0.417907,-0.314688)
(-0.299121,0.507641)
(-0.223375,0.369053)
(0.246137,-0.391641)

       */

      std::cerr << "TEST:" <<  std::endl;

      /*
      for(int kx=-3; kx<=0; kx++) {
        int lower_lim = (kx!=0)?-3:0;
        for(int ky=lower_lim; ky<3; ky++) {
          for(int kz=lower_lim; kz<3; kz++) {
            auto test_value = std::complex<T>(kx+ky+kz, 2*kx+3*ky+4*kz);
            fourier::setFourierCoefficient(linearOverdensityField, kx, ky, kz, test_value.real(), test_value.imag());
            auto check_value = fourier::getFourierCoefficient(linearOverdensityField, kx, ky, kz);
            assert(abs(test_value.real()-check_value.real())<1e-7);
            assert(abs(test_value.imag()-check_value.imag())<1e-7);
          }
        }
      }

      for(int kx=-3; kx<=0; kx++) {
        int lower_lim = (kx!=0)?-3:0;
        for(int ky=lower_lim; ky<3; ky++) {
          for(int kz=lower_lim; kz<3; kz++) {
            auto test_value = std::complex<T>(kx+ky+kz, 2*kx+3*ky+4*kz);
            auto check_value = fourier::getFourierCoefficient(linearOverdensityField, kx, ky, kz);
            assert(abs(test_value.real()-check_value.real())<1e-7);
            assert(abs(test_value.imag()-check_value.imag())<1e-7);
          }
        }
      }
       */

      // copy three times to start assembling the vx, vy, vz fields
      TField psift1k(const_cast<Grid<T> &>(grid));
      TField psift2k(const_cast<Grid<T> &>(grid));
      TField psift3k(const_cast<Grid<T> &>(grid));




      T kfft;
      size_t idx;

      const T kw = 2. * M_PI / grid.boxsize;
      const int nyquist = fourier::getNyquistModeThatMustBeReal(grid);

      fourier::applyTransformationInFourierBasis<T>(linearOverdensityField,
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
        if(abs(iix)==nyquist)
          result_x=0;
        if(abs(iiy)==nyquist)
          result_y=0;
        if(abs(iiz)==nyquist)
          result_z=0;

        return std::make_tuple(result_x, result_y, result_z);
      },
      psift1k, psift2k, psift3k);

      // TEMP DEBUG

      /*
      Field<complex<T>, T> fieldToWrite = fourier::getComplexFourierField(psift1k);
      int n = psift1k.getGrid().size;
      const int dim[3] = {n,n,n};
      numpy::SaveArrayAsNumpy("temp_debug.npy", false, 3, dim, fieldToWrite.getDataVector().data());
       */

/*
      #pragma omp parallel for schedule(static) default(shared) private(iix, iiy, iiz, kfft, idx)
      for (size_t ix = 0; ix < size; ix++) {
        pb.tick();
        for (size_t iy = 0; iy < size; iy++) {
          for (size_t iz = 0; iz < size; iz++) {

            idx = (ix * size + iy) * size + iz;

            if (ix > size / 2) iix = static_cast<int>(ix) - size; else iix = ix;
            if (iy > size / 2) iiy = static_cast<int>(iy) - size; else iiy = iy;
            if (iz > size / 2) iiz = static_cast<int>(iz) - size; else iiz = iz;

            kfft = (T) (iix * iix + iiy * iiy + iiz * iiz);

            psift1k[idx]*=T(iix)/(kfft*kw);
            psift2k[idx]*=T(iiy)/(kfft*kw);
            psift3k[idx]*=T(iiz)/(kfft*kw);
          }
        }
      }
*/
      set_zero(psift1k[0]);
      set_zero(psift2k[0]);
      set_zero(psift3k[0]);


      psift1k.toReal();
      psift2k.toReal();
      psift3k.toReal();

      // TODO: following copy operation only necessary for complex implementation

      pOff_x = std::make_shared<TRealField>(grid,false);
      pOff_y = std::make_shared<TRealField>(grid,false);
      pOff_z = std::make_shared<TRealField>(grid,false);

      auto & offx_data = pOff_x->getDataVector();
      auto & offy_data = pOff_y->getDataVector();
      auto & offz_data = pOff_z->getDataVector();

      //apply ZA:
#pragma omp parallel for schedule(static) default(shared) private(idx)
      for (size_t ix = 0; ix < size; ix++) {
        pb.tick();
        for (size_t iy = 0; iy < size; iy++) {
          for (size_t iz = 0; iz < size; iz++) {

            idx = (ix * size + iy) * size + iz;

            // position offset in physical coordinates
            offx_data[idx] = real_part_if_complex(psift1k[idx]);
            offy_data[idx] = real_part_if_complex(psift2k[idx]);
            offz_data[idx] = real_part_if_complex(psift3k[idx]);
          }
        }
      }
    }

  public:


    ZeldovichParticleGenerator(TField &linearOverdensityField,
                               const CosmologicalParameters<T> &cosmology) :
      linearOverdensityField(linearOverdensityField),
      cosmology(cosmology),
      ParticleGenerator<GridDataType>(linearOverdensityField.getGrid())
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

    void addFieldFromDifferentGridWithFilter(ZeldovichParticleGenerator &source, const Filter<T> &filter) {
      pOff_x->addFieldFromDifferentGridWithFilter(*source.pOff_x, filter);
      pOff_y->addFieldFromDifferentGridWithFilter(*source.pOff_y, filter);
      pOff_z->addFieldFromDifferentGridWithFilter(*source.pOff_z, filter);
    }

    void applyFilter(const Filter<T> &filter) {
      pOff_x->applyFilter(filter);
      pOff_y->applyFilter(filter);
      pOff_z->applyFilter(filter);
    }

    void toReal() {
      pOff_x->toReal();
      pOff_y->toReal();
      pOff_z->toReal();
    }

    virtual T getMass(const Grid<T> & onGrid) const override {
      return boxMass*onGrid.cellMassFrac;
    }

    virtual T getEps(const Grid<T> & onGrid) const override  {
      return onGrid.dx * onGrid.cellSofteningScale * 0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    auto getOffsetFields() {
      return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    auto getOffsetFields() const {
      return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    virtual particle::Particle<T> getParticleNoOffset(const Grid<T> &onGrid, size_t id) const override {
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

    virtual particle::Particle<T> getParticleNoWrap(const Grid<T> &onGrid, size_t id) const override {
      auto particle = getParticleNoOffset(onGrid, id);
      auto centroid = onGrid.getCellCentroid(id);
      particle.pos+=centroid;
      return particle;
    }



  };
}

#endif //IC_ZELDOVICH_HPP
