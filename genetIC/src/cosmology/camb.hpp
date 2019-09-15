#ifndef _CAMB_HPP_INCLUDED
#define _CAMB_HPP_INCLUDED

#include <cmath>
#include <memory>
#include <utility>
#include <map>

#include "src/cosmology/parameters.hpp"
#include "src/tools/numerics/interpolation.hpp"
#include "src/io/input.hpp"
#include "src/simulation/particles/particle.hpp"

/*!
    \namespace cosmology
    \brief Describe the cosmology adopted (parameters and transfer function)

    Allow to use cosmological parameters and to draw random fields with a cosmological power spectrum.
    The input transfer function is tied to CAMB format.
 */
namespace cosmology {

  template<typename F>
  struct CacheKeyComparator {
    bool operator()(const F &a, const F &b) const {
      bool ptr_less = std::owner_less<typename F::first_type>()(a.first, b.first);
      return ptr_less || (!ptr_less && a.second < b.second); // equivalent to normal std::pair sorting
    }
  };

  /*! \class PowerSpectrum
  * \brief Abstract base class for power spectrum calculations.
  *
  * For practical runs, the CAMB child class performs the key work. For testing runs, a 
  * PowerLawPowerSpectrum implementation is also provided.
  */
  template<typename DataType>
  class PowerSpectrum {
  public:
    using CoordinateType = tools::datatypes::strip_complex<DataType>;
    using FieldType = fields::Field<DataType, CoordinateType>;

  protected:

    using CacheKeyType = std::pair<std::weak_ptr<const grids::Grid<CoordinateType>>, particle::species>;

    //! A cache for previously calculated covariances. The key is a pair: (weak pointer to the grid, transfer fn)
    mutable std::map<CacheKeyType, std::shared_ptr<FieldType>,
      CacheKeyComparator<CacheKeyType>> calculatedCovariancesCache;

  public:

    virtual ~PowerSpectrum() {}

    //! \brief Evaluate power spectrum for a given species at wavenumber k (Mpc/h), including the normalisation
    //! transferType specifies whether to use the DM or baryon transfer function
    //! NB - for efficiency, we don't check that transferType casts to an int within bounds. This is handled by
    //! getPowerSpectrumForGrid, which does a single check before applying a guaranteed safe
    //! value to all fourier cells.
    virtual CoordinateType operator()(CoordinateType k, particle::species transferType) const = 0;

    //! Get the theoretical power spectrum appropriate for a given grid. This may be a cached copy if previously calculated.
    std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGrid(const std::shared_ptr<const grids::Grid<CoordinateType>> &grid,
                            particle::species transferType = particle::species::dm) const {

      auto cacheKey = std::make_pair(std::weak_ptr<const grids::Grid<CoordinateType>>(grid), transferType);

      auto result = this->calculatedCovariancesCache[cacheKey];

      if (result == nullptr) {
        result = getPowerSpectrumForGridUncached(grid, transferType);
        this->calculatedCovariancesCache[cacheKey] = result;
      }

      return result;
    }


  protected:
    //! Calculate the theoretical power spectrum for a given grid
    virtual std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGridUncached(std::shared_ptr<const grids::Grid<CoordinateType>> grid,
                                    particle::species transferType = particle::species::dm) const {

      CoordinateType norm = this->getPowerSpectrumNormalizationForGrid(*grid);

      auto P = std::make_shared<fields::Field<DataType, CoordinateType>>(*grid, true);

      P->forEachFourierCell([norm, this, transferType]
                              (std::complex<CoordinateType>, CoordinateType kx, CoordinateType ky,
                               CoordinateType kz) {
        CoordinateType k = sqrt(kx * kx + ky * ky + kz * kz);
        auto spec = std::complex<CoordinateType>((*this)(k, transferType) * norm, 0);
        return spec;
      });

      return P;

    }



  public:


    //! Return the box- and fft-dependent part of the normalisation of the power spectrum
    static CoordinateType getPowerSpectrumNormalizationForGrid(const grids::Grid<CoordinateType> &grid) {

      CoordinateType kw = 2. * M_PI / grid.thisGridSize;
      CoordinateType norm = kw * kw * kw / powf(2.f * M_PI, 3.f); //since kw=2pi/L, this is just 1/V_box

      // This factor Ncells was first needed when FFT normalisation changed from 1/N to 1/sqrt(N). This knowledge was previously
      // incorporated as a normalisation of the random draw rather than to the power spectrum. It makes more sense to have
      // it as a PS normalisation and restores coherence between the P(k) estimated from the field (e.g. delta dagger * delta)
      // and the theoretical P(k) calculated here. MR 2018
      CoordinateType fft_normalisation = grid.size3;

      return norm * fft_normalisation;
    }



  };

  /*! \class PowerLawPowerSpectrum
   * \brief Pure power law power spectrum, for testing purposes only.
   *
   * ns is interpreted as the power law index, and sigma8 is ignored
   */
  template<typename DataType>
  class PowerLawPowerSpectrum : public PowerSpectrum<DataType> {
  protected:
    using typename PowerSpectrum<DataType>::CoordinateType;
    CoordinateType ns;
    CoordinateType amplitude;

  public:
    PowerLawPowerSpectrum(const CosmologicalParameters<CoordinateType> &cosmology, CoordinateType amplitude)
    : ns(cosmology.ns), amplitude(amplitude) { }

    CoordinateType operator()(CoordinateType k, particle::species transferType) const override {
      if(k==0)
        return 0;
      else
        return amplitude * powf(k, ns);
    }

  };

  /*! \class CAMB
  * \brief Provides power spectra by using transfer functions from CAMB output
  */
  template<typename DataType>
  class CAMB : public PowerSpectrum<DataType> {
    using typename PowerSpectrum<DataType>::CoordinateType;

  protected:
    std::vector<CoordinateType> kInterpolationPoints; //!< Wavenumbers read from CAMB file
    std::map<particle::species, std::vector<CoordinateType>> speciesToInterpolationPoints; //!< Vector to store transfer functions
    
    const std::map<particle::species, size_t> speciesToCambColumn
      {{particle::species::dm, 1}, 
       {particle::species::baryon, 2}}; 
      //!< Columns of CAMB that we request for DM and baryons respectively
      
    std::map<particle::species, tools::numerics::LogInterpolator<CoordinateType>> speciesToTransferFunction; //!< Interpolation functions:
    CoordinateType amplitude; //!< Amplitude of the initial power spectrum
    CoordinateType ns;        //!< tensor to scalar ratio of the initial power spectrum
    mutable CoordinateType kcamb_max_in_file; //!< Maximum CAMB wavenumber. If too small compared to grid resolution, Meszaros solution will be computed

  public:
    //! Import data from CAMB file and initialise the interpolation functions used to compute the transfer functions:
    CAMB(const CosmologicalParameters<CoordinateType> &cosmology, const std::string &filename) {
      readLinesFromCambOutput(filename);
      for (auto i = speciesToInterpolationPoints.begin(); i != speciesToInterpolationPoints.end(); ++i) {
        this->speciesToTransferFunction[i->first].initialise(kInterpolationPoints, i->second);
      }
      ns = cosmology.ns;
      calculateOverallNormalization(cosmology);
    }

    CoordinateType operator()(CoordinateType k, particle::species transferType) const override {
      CoordinateType linearTransfer;
      if (k != 0)
        linearTransfer = speciesToTransferFunction.at(transferType)(k);
      else
        linearTransfer = 0.0;

      if (k > kcamb_max_in_file) {
        kcamb_max_in_file = std::numeric_limits<CoordinateType>().max();
      }

      return amplitude * powf(k, ns) * linearTransfer * linearTransfer;
    }



  protected:

      //! \brief This function imports data from a CAMB file, supplied as a file-name string argument (incamb).
      //! Both pre-2015 and post-2015 formats can be used, and the function will detect which.
      void readLinesFromCambOutput(std::string incamb) {
        kInterpolationPoints.clear();
        speciesToInterpolationPoints.clear();

        // Dealing with the update of the CAMB TF output. Both are kept for backward compatibility.
        const int c_old_camb = 7; // number of columns in camb transfer function pre 2015
        const int c_new_camb = 13; // number of columns in camb transfer function post 2015
        int numCols;
        size_t j;

        // Import data from CAMB file:
        std::vector<double> input;

        // Have to do this while stripping out any column headers that might be present in the file:
        io::getBufferIgnoringColumnHeaders(input, incamb);

        // Check whether the input file is in the pre-2015 or post-2015 format (and throw an error if it is neither).
        numCols = io::getNumberOfColumns(incamb);

        if (numCols == c_old_camb) {
          std::cerr << "Using pre 2015 CAMB transfer function" << std::endl;
        } else if (numCols == c_new_camb) {
          std::cerr << "Using post 2015 CAMB transfer function" << std::endl;
        } else {
          throw std::runtime_error("CAMB transfer file doesn't have a sensible number of rows and columns");
        }


        CoordinateType transferNormalisation = input[1]; // to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural
        // Copy file into vectors. Normalise both so that the DM tranfer function starts at 1.
        for (j = 0; j < input.size() / numCols; j++) {
          if (input[numCols * j] > 0) {
            // hard-coded to first two columns of CAMB file -
            kInterpolationPoints.push_back(CoordinateType(input[numCols * j]));
            for (auto i = speciesToCambColumn.begin(); i != speciesToCambColumn.end(); ++i) {
              speciesToInterpolationPoints[i->first].push_back(CoordinateType(input[numCols * j + i->second]) / transferNormalisation);
            }
          } else continue;
        }

        // Extend high-k range using Meszaros solution
        // This is a very naive approximation and a big warning will be issued if the power is actually evaluated
        // at these high k's (see operator() below).


        kcamb_max_in_file = kInterpolationPoints.back();
        CoordinateType keq = 0.01;
        while (kInterpolationPoints.back() < 1e7) {
          kInterpolationPoints.push_back(kInterpolationPoints.back() * 1.1);
          CoordinateType kratio = kInterpolationPoints.back() / kcamb_max_in_file;
          for (auto i = speciesToInterpolationPoints.begin(); i != speciesToInterpolationPoints.end(); ++i) {
            i->second.push_back(i->second.back() * pow(kratio, -2.0) * log(kInterpolationPoints.back() / keq) / log(kcamb_max_in_file / keq));
          }
        }



      }

    //! Calculate the theoretical power spectrum for a given grid
    std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGridUncached(std::shared_ptr<const grids::Grid<CoordinateType>> grid,
                                    particle::species transferType = particle::species::dm) const override {

      auto P = PowerSpectrum<DataType>::getPowerSpectrumForGridUncached(grid, transferType);

      if (kcamb_max_in_file == std::numeric_limits<CoordinateType>().max()) {
        std::cerr << "WARNING: maximum k in CAMB input file is insufficient" << std::endl;
        std::cerr << "         extrapolating using naive Meszaros solution" << std::endl;
      }

      return P;

    }

      //! Return the cosmology-dependent part of the normalisation of the power spectrum.
      void calculateOverallNormalization(const CosmologicalParameters<CoordinateType> &cosmology) {
        CoordinateType ourGrowthFactor = growthFactor(cosmology);
        CoordinateType growthFactorNormalized = ourGrowthFactor / growthFactor(cosmologyAtRedshift(cosmology, 0));
        CoordinateType sigma8PreNormalization = calculateLinearVarianceInSphere(8.);
        CoordinateType linearRenormFactor = (cosmology.sigma8 / sigma8PreNormalization) * growthFactorNormalized;

        amplitude = linearRenormFactor * linearRenormFactor;
      }


      //! Compute the variance in a spherical top hat window
      CoordinateType calculateLinearVarianceInSphere(CoordinateType radius,
                                                     particle::species transferType = particle::species::dm) const {

        CoordinateType s = 0., k, t;

        CoordinateType amp = 9. / 2. / M_PI / M_PI;
        CoordinateType kmax = std::min(kInterpolationPoints.back(), 200.0 / radius)*0.999999;
        CoordinateType kmin = kInterpolationPoints[0]*1.000001;

        CoordinateType dk = (kmax - kmin) / 50000.;
        auto &interpolator = this->speciesToTransferFunction.at(transferType);
        for (k = kmin; k < kmax; k += dk) {


          t = interpolator(k);

          // Multiply power spectrum by the fourier transform of the spherical top hat, to give the fourier transform
          // of the averaged (convolved) power spectrum over the sphere.
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
