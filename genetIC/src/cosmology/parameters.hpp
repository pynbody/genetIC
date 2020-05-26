#ifndef _COSMO_HPP_INCLUDED
#define _COSMO_HPP_INCLUDED

#include <functional>
#include <climits>


namespace cosmology {

  template<typename DataType>
  class CAMB;

  /*! \struct CosmologicalParameters
      \brief Stores data about the cosmological model being assumed.

      This allow us to quickly pass a long list of cosmological parameters to
      different parts of the code.
  */
  template<typename FloatType>
  struct CosmologicalParameters {
    FloatType OmegaM0; //!< Total matter density fraction (dark matter + baryons);
    FloatType OmegaLambda0; //!< Dark energy density fraction
    FloatType OmegaBaryons0; //!< Baryon density fraction
    FloatType hubble; //!< Hubble rate in units of 100kms^{-1}Mpc^{-1}
    FloatType redshift; //!< Redshift at which we are generating the initial conditions.
    FloatType scalefactor; //!< Scale factor at this redshift, relative to the scale factor at z = 0
    FloatType scalefactorAtDecoupling; //!< Scale factor at baryon-photon decoupling (z~150, NOT photon-baryon decoupling at z~1100)
    FloatType sigma8; //!< Sigma8 parameter
    FloatType ns; //!< Scalar spectral index
    FloatType TCMB; //!< CMB temperature today, in K
  };

  //! Computes an estimate of the linear growth factor.
  template<typename FloatType>
  FloatType growthFactor(const CosmologicalParameters<FloatType> &cosmology) {
    const FloatType a = cosmology.scalefactor;
    const FloatType Om = cosmology.OmegaM0;
    const FloatType Ol = cosmology.OmegaLambda0;

    // Square of the Hubble rate as a fraction of the present day, Hubble rate, ie, (H/H0)^2:
    // NB - this uses the curvature term, which isn't (currently) used elsewhere in the code.
    FloatType Hsq = cosmology.OmegaM0 / powf(a, 3.0) + (1. - Om - Ol) / a / a + Ol;
    // Structure growth factor, analytic approximation (exact requires solding an ode):
    FloatType d = 2.5 * a * Om / pow(a, 3.0) / Hsq / (pow(Om / Hsq / a / a / a, 4.f / 7.f) - Ol / Hsq +
                                                       (1. + 0.5 * Om / pow(a, 3.0) / Hsq) *
                                                       (1. + 1. / 70. * Ol / Hsq));

    return d;
  }


  //! Returns a copy of the cosmological parameters with the scalefactor updated to match redshift z
  template<typename FloatType>
  CosmologicalParameters<FloatType>
  cosmologyAtRedshift(const CosmologicalParameters<FloatType> &referenceCosmology, float redshift) {
    CosmologicalParameters<FloatType> retVal(referenceCosmology);
    retVal.scalefactor = 1. / (1. + redshift);
    return retVal;
  }

  //! Calculate the ratio between the velocity and the position offset in internal units
  template<typename FloatType>
  FloatType zeldovichVelocityToOffsetRatio(const CosmologicalParameters<FloatType> &cosmology) {

    // TODO: hardcoded value of f=1 is inaccurate - should be a function of omega
    FloatType f = 1.0;
    // According to Carrol and Press (1992) should use:
    // f = pow(( Om/(a*a*a) )/( Om/(a*a*a) + (1.0 - Om - Ol)/(a*a) + Ol ),4.0/7.0);

    // For ease of reading the expression below:
    FloatType Om = cosmology.OmegaM0;
    FloatType Ol = cosmology.OmegaLambda0;
    FloatType a = cosmology.scalefactor;


    FloatType velocityToOffsetRatio = f * 100. * sqrt(Om / (a * a * a) + Ol) * sqrt(a);
    // this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)

    // Could also include curvature term, ie, f * 100. * sqrt( Om/(a*a*a) + (1 - Om - Ol)/(a*a) + Ol)*sqrt(a);
    // Seems a bit inconsistent that this is used in growthFactor, implicitly, in the definition of Hsq

    return velocityToOffsetRatio;


  }


  //! \brief Dump an estimated power spectrum for the field, alongside the specified theory power spectrum, to disk
  // TODO Refactor this to use grid methods and normalisations method. It could be compacted in a few lines method
  // and avoid repeating assumptions from elsewhere in the code.
  template<typename DataType, typename FloatType=tools::datatypes::strip_complex<DataType>>
  void dumpPowerSpectrum(const fields::Field<DataType> &field,
                         const fields::Field<DataType> &P0, const std::string &filename) {

    // Strategy here is to estimate the power spectrum by summing over all points in the
    // generated field in Fourier space, and assigning them to fixed-width k bins according
    // to k = \sqrt{kx^2 + ky^2 + kz^2}. Take the average of |\delta_k|^2 in each bin to get
    // the power spectrum at the average value of k in that bin. More accurate for high k-modes, since
    // there are more of them able to fit in the simulation box.

    field.ensureFourierModesAreMirrored();
    P0.ensureFourierModesAreMirrored();

    int res = field.getGrid().size; // Over-density field
    int nBins = 100; // Bins used to estimate the power spectrum
    std::vector<FloatType> inBin(nBins); // Wavenumbers contributing to a given k bin.
    std::vector<FloatType> kbin(nBins); // Sum of k = \sqrt{kx^2 + ky^2 + kz^2} value contributing to this bin.
    std::vector<FloatType> Gx(nBins); // Sum of |\delta_k|^2 contributing to this bin.
    std::vector<FloatType> Px(nBins); // Sum of 'exact' values of power spectrum for k contributing to this bin.


    // Get the range of k-modes to compute the power spectrum for, based on the size of the simulation box
    // and grid-scale:
    const FloatType boxLength = field.getGrid().thisGridSize;
    FloatType kmax = M_PI / boxLength * (FloatType) res, kmin = 2.0f * M_PI / boxLength, dklog =
      log10(kmax / kmin) / nBins, kw = 2.0f * M_PI / boxLength;

    // Initialise storage for bins:
    int ix, iy, iz, idx;
    FloatType kfft;

    for (ix = 0; ix < nBins; ix++) {
      inBin[ix] = 0;
      Gx[ix] = 0.0f;
      kbin[ix] = 0.0f;
    }


    // Iterate over all Fourier modes of the field:
    for (ix = -res / 2; ix < res / 2 + 1; ix++)
      for (iy = -res / 2; iy < res / 2 + 1; iy++)
        for (iz = -res / 2; iz < res / 2 + 1; iz++) {
          // Compute square of the Fourier mode:
          auto fieldValue = field.getFourierCoefficient(ix, iy, iz);
          FloatType vabs = std::abs(fieldValue);
          vabs *= vabs;

          // Compute k for this Fourier mode:
          kfft = sqrt(ix * ix + iy * iy + iz * iz);
          FloatType k = kfft * kw;

          // .. logarithmic spacing in k
          idx = (int) ((1.0f / dklog * log10(k / kmin)));

          // Make sure we don't go outside the bins by mistake:
          if (k >= kmin && k < kmax) {

            Gx[idx] += vabs; // Sum squares of the field
            Px[idx] += P0.getFourierCoefficient(ix, iy, iz).real(); // Sum 'exact' values of power spectrum in this bin.
            kbin[idx] += k; // Sum of k contributing to this bin.
            inBin[idx]++; // Total number in this bin

          } else { continue; }

        }


    // ... convert to comoving units ...
    std::ofstream ofs(filename);
    FloatType psnorm =
      1 / (CAMB<FloatType>::getPowerSpectrumNormalizationForGrid(field.getGrid()));
    // Pre-v1 this had an extra factor pow((2.0 * M_PI), 3.0), which made it differ from CAMB normalisation
    // Updated to prevent confusion.

    for (ix = 0; ix < nBins; ix++) {

      if (inBin[ix] > 0) {


        ofs << std::setw(16) << pow(10., log10(kmin) + dklog * (ix + 0.5)) // Middle of the k-bin
            << std::setw(16) << kbin[ix] / inBin[ix] // Average value of k in this bin
            << std::setw(16)
            << (FloatType) (Px[ix] / inBin[ix]) * psnorm // Average of exact power spectrum for k in this bin
            << std::setw(16) << (FloatType) (Gx[ix] / inBin[ix]) *
                                psnorm // Average of |delta_k|^2 in this bin, ie, estimated power spectrum.
            << std::setw(16) << inBin[ix] // Number in this bin
            << std::endl;

      }
    }
    ofs.close();


  }

  //! Convert a density field to a potential field, in-place.
  template<typename DataType, typename FloatType=tools::datatypes::strip_complex<DataType>>
  void densityToPotential(fields::Field<DataType, FloatType> &field, const CosmologicalParameters<FloatType> &cosmo)
  {
    field.toFourier();

    const FloatType prefac =
      3. / 2. * cosmo.OmegaM0 / cosmo.scalefactor
      * 100. * 100. / (299792.) /
      (299792.); // =3/2 Om0/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)

    field.forEachFourierCell([&prefac](complex<FloatType> val, FloatType k1, FloatType k2, FloatType k3) {
      FloatType ksq = k1 * k1 + k2 * k2 + k3 * k3;
      FloatType ksq_inv = 1.0 / ksq;
      if(k2==0.0) ksq_inv = 0.0;
      return -prefac * val * ksq_inv;
    });

  }

}

#endif
