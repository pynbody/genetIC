#ifndef _CAMB_HPP_INCLUDED
#define _CAMB_HPP_INCLUDED

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

  template<typename T>
  inline constexpr bool operator<(const std::weak_ptr<T> &a, const std::weak_ptr<T> &b) {
    return std::owner_less<std::weak_ptr<T>>()(a, b);
  }

  /*! \class CAMB
  * \brief Load in and provide interpolation routines for the cosmological power spectrum.
  *
  * Currently tied to the CAMB transfer function output format, though this could easily
  * be relaxed in future by creating an abstract base class and deriving different
  * classes for alternative approaches. */
  template<typename DataType>
  class CAMB {
  protected:
    using CoordinateType = tools::datatypes::strip_complex<DataType>;
    using FieldType = fields::Field<DataType, CoordinateType>;
    std::vector<CoordinateType> kcamb; //!< Wavenumbers read from CAMB file
    std::vector<std::vector<CoordinateType>> T; //!< Vector to store transfer functions

    /*! \brief Map from the species (baryon/dm) to the correct transfer function.

       Elements are all zero if dmOnly = true, and (dm->0, baryon->1) if dmOnly = false.
       */
    std::map<particle::species, size_t> speciesToTransferId;

    using CacheKeyType = std::pair<std::weak_ptr<const grids::Grid<CoordinateType>>, size_t>;


    //! A cache for previously calculated covariances. The key is a pair: (weak pointer to the grid, transfer fn)
    mutable std::map<CacheKeyType, std::shared_ptr<FieldType>,
      CacheKeyComparator<CacheKeyType>> calculatedCovariancesCache;

    size_t cambCols[2] = {1, 2}; //!< Columns of CAMB that we request:
    std::vector<tools::numerics::Interpolator<CoordinateType>> allInterpolators; //!< Interpolation functions:
    CoordinateType amplitude; //!< Amplitude of the initial power spectrum
    CoordinateType ns;        //!< tensor to scalar ratio of the initial power spectrum
    mutable CoordinateType kcamb_max_in_file; //!< Maximum CAMB wavenumber. If too small compared to grid resolution, Meszaros solution will be computed

  public:


    const size_t nTransfers = 2; //!< Total number of transfer functions stored:
    bool dataRead; //!< True if we have read in data.

    CAMB(bool useOnlyDM = true) {
      this->dataRead = false;
      // Initialise the interpolators and transfer functions vectors to be the
      // appropriate size, containing default constructed objects.
      if (useOnlyDM) {
        this->allInterpolators.resize(1);
        this->T.resize(1);
        speciesToTransferId[particle::species::dm] = 0;
        speciesToTransferId[particle::species::baryon] = 0;
      } else {
        this->allInterpolators.resize(this->nTransfers);
        this->T.resize(this->nTransfers);
        speciesToTransferId[particle::species::dm] = 0;
        speciesToTransferId[particle::species::baryon] = 1;
      }
    }

    //! Enables the use of the baryon transfer function
    void enableAllTransfers() {
      if (this->T.size() < 2) {
        this->kcamb.clear();
        this->T.clear();
      }
      this->T.resize(this->nTransfers);
      this->allInterpolators.resize(this->nTransfers);
      speciesToTransferId[particle::species::dm] = 0;
      speciesToTransferId[particle::species::baryon] = 1;
    }

    //! Checks whether the CAMB object actually stores any data:
    bool isUsable() const {
      return (kcamb.size() > 0);
    }

    //! Import data from CAMB file and initialise the interpolation functions used to compute the transfer functions:
    void read(const std::string &filename, const CosmologicalParameters<CoordinateType> &cosmology) {
      readLinesFromCambOutput(filename);
      for (size_t i = 0; i < this->T.size(); i++) {
        this->allInterpolators[i].initialise(kcamb, this->T[i]);
      }

      // a bit awkward that we have to copy this value:
      ns = cosmology.ns;

      calculateOverallNormalization(cosmology);

      this->dataRead = true;
    }

  protected:
    //! \brief This function imports data from a CAMB file, supplied as a file-name string argument (incamb).
    //! Both pre-2015 and post-2015 formats can be used, and the function will detect which.
    void readLinesFromCambOutput(std::string incamb) {
      kcamb.clear();
      for (size_t i = 0; i < this->T.size(); i++) {
        T[i].clear();
      }


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
          kcamb.push_back(CoordinateType(input[numCols * j]));
          for (size_t i = 0; i < this->T.size(); i++) {
            T[i].push_back(CoordinateType(input[numCols * j + this->cambCols[i]]) / transferNormalisation);
          }
        } else continue;
      }

      // Extend high-k range using Meszaros solution
      // This is a very naive approximation and a big warning will be issued if the power is actually evaluated
      // at these high k's (see operator() below).

      std::vector<CoordinateType> T_f;
      for (size_t i = 0; i < this->T.size(); i++) {
        T_f.push_back(T[i].back());
      }
      kcamb_max_in_file = kcamb.back();
      CoordinateType keq = 0.01;
      while (kcamb.back() < 1e7) {
        kcamb.push_back(kcamb.back() + 1.0);
        CoordinateType kratio = kcamb.back() / kcamb_max_in_file;
        for (size_t i = 0; i < this->T.size(); i++) {
          T[i].push_back(T_f[i] * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
        }
      }


    }

  public:

    //! \brief Evaluate power spectrum nTransfer at wavenumber k (Mpc/h), including the normalisation
    //! transferType specifies whether to use the DM or baryon transfer function
    //! NB - for efficiency, we don't check that transferType casts to an int within bounds. This is handled by
    //! getPowerSpectrumForGrid, which does a single check before applying a guaranteed safe
    //! value to all fourier cells.
    CoordinateType operator()(CoordinateType k, particle::species transferType) const {
      CoordinateType linearTransfer;
      if (k != 0)
        linearTransfer = allInterpolators[speciesToTransferId.at(transferType)](k);
        // speciesToTransferId ensures that if extra transfer functions are switched off, then we always
        // call the dark matter allInterpolators, rather than falling out of bounds.
      else
        linearTransfer = 0.0;

      if (k > kcamb_max_in_file) {
        kcamb_max_in_file = std::numeric_limits<CoordinateType>().max();
      }

      return amplitude * powf(k, ns) * linearTransfer * linearTransfer;
    }

    //! Get the theoretical power spectrum appropriate for a given grid. This may be a cached copy if previously calculated.
    std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGrid(const std::shared_ptr<const grids::Grid<CoordinateType>> &grid,
                            particle::species transferType = particle::species::dm) const {
      size_t transferId = speciesToTransferId.at(transferType);

      auto cacheKey = std::make_pair(std::weak_ptr<const grids::Grid<CoordinateType>>(grid), transferId);

      auto result = this->calculatedCovariancesCache[cacheKey];

      if (result == nullptr) {
        result = getPowerSpectrumForGridUncached(grid, transferType);
        this->calculatedCovariancesCache[cacheKey] = result;
      }

      return result;

    }


  protected:
    //! Calculate the theoretical power spectrum for a given grid
    std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGridUncached(std::shared_ptr<const grids::Grid<CoordinateType>> grid,
                                    particle::species transferType = particle::species::dm) const {
      assert(kcamb.size() == T[0].size());

      CoordinateType norm = getPowerSpectrumNormalizationForGrid(*grid);

      auto P = std::make_shared<fields::Field<DataType, CoordinateType>>(*grid, true);

      P->forEachFourierCell([norm, this, transferType]
                              (std::complex<CoordinateType>, CoordinateType kx, CoordinateType ky,
                               CoordinateType kz) {
        CoordinateType k = sqrt(kx * kx + ky * ky + kz * kz);
        auto spec = std::complex<CoordinateType>((*this)(k, transferType) * norm, 0);
        return spec;
      });

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


    //! Compute the variance in a spherical top hat window
    CoordinateType calculateLinearVarianceInSphere(CoordinateType radius,
                                                   particle::species transferType = particle::species::dm) const {

      CoordinateType s = 0., k, t;

      CoordinateType amp = 9. / 2. / M_PI / M_PI;
      CoordinateType kmax = std::min(kcamb.back(), 200.0 / radius);
      CoordinateType kmin = kcamb[0];

      CoordinateType dk = (kmax - kmin) / 50000.;
      auto &interpolator = this->allInterpolators[speciesToTransferId.at(transferType)];
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
