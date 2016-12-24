//
// Created by Andrew Pontzen on 20/12/2016.
//

#ifndef IC_ZELDOVICH_HPP
#define IC_ZELDOVICH_HPP

#include <complex>
#include "../cosmo.hpp"
#include "../progress/progress.hpp"
#include "../field.hpp"
#include "generator.hpp"

template<typename T>
struct CosmologicalParameters;

namespace particle {

  template<typename T>
  class ZeldovichParticleGenerator : public ParticleGenerator<T> {
  protected:
    using TField = std::vector<std::complex<T>>;
    using TRealField = Field<T,T>;

    T velocityToOffsetRatio, particleMass;
    Field<complex<T>, T> &linearOverdensityField;
    const CosmologicalParameters<T> &cosmology;
    using ParticleGenerator<T>::grid;

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
    void calculateParticleMass() {
      particleMass = 27.78 * cosmology.OmegaM0 * powf(grid.boxsize / (T) (grid.size), 3.0);
    }
    virtual void calculateOffsetFields() {
      // TODO: refactorise this horrible long method

      size_t size = grid.size;
      size_t size3 = grid.size3;
      progress::ProgressBar pb("zeldovich", size * 2);

      // make three arrays for manipulating in fourier space
      auto psift1k = std::vector<std::complex<T>>(size3, 0);
      auto psift2k = std::vector<std::complex<T>>(size3, 0);
      auto psift3k = std::vector<std::complex<T>>(size3, 0);

      // get a reference to the density field in fourier space
      linearOverdensityField.toFourier();
      auto &pField_k = linearOverdensityField.getDataVector();

      int iix, iiy, iiz;
      T kfft;
      size_t idx;

      T kw = 2. * M_PI / grid.boxsize;

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

            psift1k[idx].real(-pField_k[idx].imag() / (T) (kfft) * iix / kw);
            psift1k[idx].imag(pField_k[idx].real() / (T) (kfft) * iix / kw);
            psift2k[idx].real(-pField_k[idx].imag() / (T) (kfft) * iiy / kw);
            psift2k[idx].imag(pField_k[idx].real() / (T) (kfft) * iiy / kw);
            psift3k[idx].real(-pField_k[idx].imag() / (T) (kfft) * iiz / kw);
            psift3k[idx].imag(pField_k[idx].real() / (T) (kfft) * iiz / kw);
          }
        }
      }

      psift1k[0] = complex<T>(0., 0.);
      psift2k[0] = complex<T>(0., 0.);
      psift3k[0] = complex<T>(0., 0.);

      fft(psift1k.data(), psift1k.data(), size,
          -1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
      fft(psift2k.data(), psift2k.data(), size, -1); //same
      fft(psift3k.data(), psift3k.data(), size, -1); //same

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
            offx_data[idx] = psift1k[idx].real();
            offy_data[idx] = psift2k[idx].real();
            offz_data[idx] = psift3k[idx].real();
          }
        }
      }
    }

  public:


    ZeldovichParticleGenerator(Field<complex<T> ,T> &linearOverdensityField,
                               const CosmologicalParameters<T> &cosmology) :
      linearOverdensityField(linearOverdensityField),
      cosmology(cosmology),
      ParticleGenerator<T>(linearOverdensityField.getGrid())
    {
      recalculate();
    }

    void recalculate() override {
      calculateVelocityToOffsetRatio();
      calculateParticleMass();
      calculateOffsetFields();
    }

    void addFieldFromDifferentGrid(const ZeldovichParticleGenerator &source) {
      pOff_x->addFieldFromDifferentGrid(*source.pOff_x);
      pOff_y->addFieldFromDifferentGrid(*source.pOff_y);
      pOff_z->addFieldFromDifferentGrid(*source.pOff_z);
    }

    virtual T getMass() const override {
      return particleMass;
    }

    virtual T getEps() const override  {
      return grid.dx * 0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    auto getOffsetFields() {
      return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    auto getOffsetFields() const {
      return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    virtual particle::Particle<T> getParticleNoOffset(size_t id) const override {
      particle::Particle<T> particle;

      particle.pos.x = (*pOff_x)[id];
      particle.pos.y = (*pOff_y)[id];
      particle.pos.z = (*pOff_z)[id];

      particle.vel.x = (*pOff_x)[id] * velocityToOffsetRatio;
      particle.vel.y = (*pOff_y)[id] * velocityToOffsetRatio;
      particle.vel.z = (*pOff_z)[id] * velocityToOffsetRatio;

      particle.mass = getMass();
      particle.soft = getEps();

      return particle;
    }

    virtual void getParticleFromOffset(particle::Particle<T> &particle) const override {

      particle.vel.x = grid.getFieldInterpolated(particle.pos, *pOff_x);
      particle.vel.y = grid.getFieldInterpolated(particle.pos, *pOff_y);
      particle.vel.z = grid.getFieldInterpolated(particle.pos, *pOff_z);

      particle.pos+=particle.vel;
      particle.vel*=velocityToOffsetRatio;

      grid.simWrap(particle.pos);

      particle.mass = getMass();
      particle.soft = getEps();
    }

    virtual particle::Particle<T> getParticleNoWrap(size_t id) const override {
      auto particle = getParticleNoOffset(id);
      auto centroid = grid.getCellCentroid(id);
      particle.pos+=centroid;
      return particle;
    }


  };
}

#endif //IC_ZELDOVICH_HPP
