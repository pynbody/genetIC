#ifndef _CAMB_HPP_INCLUDED
#define _CAMB_HPP_INCLUDED

#include "parameters.hpp"
#include "../numerics/interpolation.hpp"

namespace cosmology {

  /** Load in and provide interpolation routines for the cosmological power spectrum.
   *
   * Currently tied to the CAMB transfer function output format, though this could easily
   * be relaxed in future by creating an abstract base class and deriving different
   * classes for alternative approaches. */
  template<typename FloatType>
  class CAMB {
  protected:
    std::vector <FloatType> kcamb;
    std::vector <FloatType> Tcamb;
    numerics::Interpolator <FloatType> interpolator;
    FloatType amplitude;
    FloatType ns;
    mutable FloatType kcamb_max_in_file;

  public:

    bool isUsable() const {
      return (kcamb.size() > 0);
    }

    void read(std::string incamb, const CosmologicalParameters <FloatType> &cosmology) {

      readLinesFromCambOutput(incamb);
      interpolator.initialise(kcamb, Tcamb);

      // a bit awkward that we have to copy this value:
      ns = cosmology.ns;

      calculateOverallNormalization(cosmology);


      // TODO: here is where we'd insert a conversion of the power spectrum to match the real-space correlation function.
      // This work was paused because it's not so obvious whether we really want it or not (literature generally claims
      // yes, but claims are not totally convincing.)

      /*
      realspace::RealSpaceGenerators<FloatType> obj(200.0,8192*2);

      cerr << "k=" << obj.generateKArray() << endl;
      cerr << "Pk=" << obj.generatePkArray(*this) << endl;
      */

    }

  protected:
    void readLinesFromCambOutput(std::string incamb) {
      kcamb.clear();
      Tcamb.clear();

      const int c = 7; // number of columns in transfer function
      size_t j;

      std::vector<double> input;

      io::getBuffer(input, incamb);

      if (input.size() < c || input.size() % c != 0) {
        throw std::runtime_error("CAMB transfer file doesn't have a sensible number of rows and columns");
      }


      FloatType ap = input[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

      for (j = 0; j < input.size() / c; j++) {
        if (input[c * j] > 0) {
          // hard-coded to first two columns of CAMB file -
          kcamb.push_back(FloatType(input[c * j]));
          Tcamb.push_back(FloatType(input[c * j + 1]) / ap);
        } else continue;
      }

      size_t nCambLines = kcamb.size();

      // Extend high-k range using Meszaros solution
      // This is a very naive approximation and a big warning will be issued if the power is actually evaluated
      // at these high k's (see operator() below).

      FloatType Tcamb_f = Tcamb.back();
      kcamb_max_in_file = kcamb.back();
      FloatType keq = 0.01;
      while (kcamb.back() < 1000) {
        kcamb.push_back(kcamb.back() + 1.0);
        FloatType kratio = kcamb.back() / kcamb_max_in_file;
        Tcamb.push_back(Tcamb_f * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
      }


    }

  public:

    FloatType operator()(FloatType k) const {
      /* Evaluate the power spectrum at wavenumber k (Mpc/h), including the normalisation */
      FloatType linearTransfer;
      if (k != 0)
        linearTransfer = interpolator(k);
      else
        linearTransfer = 0.0;

      if (k > kcamb_max_in_file) {
        std::cerr << "WARNING: maximum k in CAMB input file is insufficient (" << kcamb_max_in_file << ")" << std::endl;
        std::cerr << "         extrapolating using naive Meszaros solution" << std::endl;
        kcamb_max_in_file = std::numeric_limits<FloatType>().max();
      }

      return amplitude * powf(k, ns) * linearTransfer * linearTransfer;
    }

    std::vector <FloatType> getPowerSpectrumForGrid(const grids::Grid<FloatType> &grid) const {
      /* Get the variance for each Fourier cell of the specified grid  */
      assert(kcamb.size() == Tcamb.size());

      FloatType norm = getPowerSpectrumNormalizationForGrid(grid);

      std::vector <FloatType> P(grid.size3);

#pragma omp parallel for
      for (size_t i = 0; i < grid.size3; ++i) {
        FloatType k = grid.getFourierCellAbsK(i);
        P[i] = (*this)(k) * norm;
      }

      return P;

    }


  protected:
    void calculateOverallNormalization(const CosmologicalParameters <FloatType> &cosmology) {
      FloatType ourGrowthFactor = growthFactor(cosmology);
      FloatType growthFactorNormalized = ourGrowthFactor / growthFactor(cosmologyAtRedshift(cosmology, 0));
      FloatType sigma8PreNormalization = calculateLinearVarianceInSphere(8.);
      FloatType linearRenormFactor = (cosmology.sigma8 / sigma8PreNormalization) * growthFactorNormalized;

      amplitude = linearRenormFactor * linearRenormFactor;

    }

    FloatType getPowerSpectrumNormalizationForGrid(const grids::Grid<FloatType> &grid) const {

      FloatType kw = 2. * M_PI / grid.boxsize;
      FloatType norm = kw * kw * kw / powf(2. * M_PI, 3.); //since kw=2pi/L, this is just 1/V_box

      return norm;
    }

  public:


    FloatType calculateLinearVarianceInSphere(FloatType radius) const {

      FloatType s = 0., k, t;

      FloatType amp = 9. / 2. / M_PI / M_PI;
      FloatType kmax = std::min(kcamb.back(), 200.0 / radius);
      FloatType kmin = kcamb[0];

      FloatType dk = (kmax - kmin) / 50000.;
      for (k = kmin; k < kmax; k += dk) {

        t = interpolator(k);

        s += powf(k, ns + 2.) *
             ((sin(k * radius) - k * radius * cos(k * radius)) / ((k * radius) * (k * radius) * (k * radius))) *
             ((sin(k * radius) - k * radius * cos(k * radius)) / ((k * radius) * (k * radius) * (k * radius))) * t * t;

      }


      s = sqrt(s * amp * dk);
      return s;

    }

  };

}

#endif