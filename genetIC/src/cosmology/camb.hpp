#ifndef _CAMB_HPP_INCLUDED
#define _CAMB_HPP_INCLUDED

#include "src/cosmology/parameters.hpp"
#include "src/tools/numerics/interpolation.hpp"

/*!
    \namespace cosmology
    \brief Describe the cosmology adopted (parameters and transfer function)

    Allow to use cosmological parameters and to draw random fields with a cosmological power spectrum.
    The input transfer function is tied to CAMB format.
 */
namespace cosmology {

  /** Load in and provide interpolation routines for the cosmological power spectrum.
   *
   * Currently tied to the CAMB transfer function output format, though this could easily
   * be relaxed in future by creating an abstract base class and deriving different
   * classes for alternative approaches. */
  template<typename DataType>
  class CAMB {
  protected:
    using CoordinateType=tools::datatypes::strip_complex<DataType>;
    std::vector<CoordinateType> kcamb;
    std::vector<CoordinateType> Tcamb;
    tools::numerics::Interpolator<CoordinateType> interpolator;
    CoordinateType amplitude;
    CoordinateType ns;
    mutable CoordinateType kcamb_max_in_file;

  public:

    bool isUsable() const {
      return (kcamb.size() > 0);
    }

    void read(std::string incamb, const CosmologicalParameters <CoordinateType> &cosmology) {

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


      CoordinateType ap = input[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

      for (j = 0; j < input.size() / c; j++) {
        if (input[c * j] > 0) {
          // hard-coded to first two columns of CAMB file -
          kcamb.push_back(CoordinateType(input[c * j]));
          Tcamb.push_back(CoordinateType(input[c * j + 1]) / ap);
        } else continue;
      }

      // Extend high-k range using Meszaros solution
      // This is a very naive approximation and a big warning will be issued if the power is actually evaluated
      // at these high k's (see operator() below).

      CoordinateType Tcamb_f = Tcamb.back();
      kcamb_max_in_file = kcamb.back();
      CoordinateType keq = 0.01;
      while (kcamb.back() < 1e7) {
        kcamb.push_back(kcamb.back() + 1.0);
        CoordinateType kratio = kcamb.back() / kcamb_max_in_file;
        Tcamb.push_back(Tcamb_f * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
      }


    }

  public:

    CoordinateType operator()(CoordinateType k) const {
      /* Evaluate the power spectrum at wavenumber k (Mpc/h), including the normalisation */
      CoordinateType linearTransfer;
      if (k != 0)
        linearTransfer = interpolator(k);
      else
        linearTransfer = 0.0;

      if (k > kcamb_max_in_file) {
        kcamb_max_in_file = std::numeric_limits<CoordinateType>().max();
      }

      return amplitude * powf(k, ns) * linearTransfer * linearTransfer;
    }

    std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGrid(const grids::Grid<CoordinateType> &grid) const {
      /* Get the variance for each Fourier cell of the specified grid  */
      assert(kcamb.size() == Tcamb.size());

      CoordinateType norm = getPowerSpectrumNormalizationForGrid(grid);

      auto P = std::make_shared<fields::Field<DataType, CoordinateType>>(grid, true);

      P->forEachFourierCell([norm, this]
                                (std::complex<CoordinateType>, CoordinateType kx, CoordinateType ky, CoordinateType kz) {
        CoordinateType k = sqrt(kx*kx+ky*ky+kz*kz);
        return std::complex<CoordinateType>((*this)(k)*norm,0);
      });

      if (kcamb_max_in_file == std::numeric_limits<CoordinateType>().max()) {
        std::cerr << "WARNING: maximum k in CAMB input file is insufficient" << std::endl;
        std::cerr << "         extrapolating using naive Meszaros solution" << std::endl;
      }

      return P;

    }


  protected:
    void calculateOverallNormalization(const CosmologicalParameters <CoordinateType> &cosmology) {
      CoordinateType ourGrowthFactor = growthFactor(cosmology);
      CoordinateType growthFactorNormalized = ourGrowthFactor / growthFactor(cosmologyAtRedshift(cosmology, 0));
      CoordinateType sigma8PreNormalization = calculateLinearVarianceInSphere(8.);
      CoordinateType linearRenormFactor = (cosmology.sigma8 / sigma8PreNormalization) * growthFactorNormalized;

      amplitude = linearRenormFactor * linearRenormFactor;

    }

    CoordinateType getPowerSpectrumNormalizationForGrid(const grids::Grid<CoordinateType> &grid) const {

      CoordinateType kw = 2. * M_PI / grid.boxsize;
      CoordinateType norm = kw * kw * kw / powf(2. * M_PI, 3.); //since kw=2pi/L, this is just 1/V_box

      return norm;
    }

  public:


    CoordinateType calculateLinearVarianceInSphere(CoordinateType radius) const {

      CoordinateType s = 0., k, t;

      CoordinateType amp = 9. / 2. / M_PI / M_PI;
      CoordinateType kmax = std::min(kcamb.back(), 200.0 / radius);
      CoordinateType kmin = kcamb[0];

      CoordinateType dk = (kmax - kmin) / 50000.;
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
