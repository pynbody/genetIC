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

   /*! \class CAMB
   * \brief Load in and provide interpolation routines for the cosmological power spectrum.
   *
   * Currently tied to the CAMB transfer function output format, though this could easily
   * be relaxed in future by creating an abstract base class and deriving different
   * classes for alternative approaches. */
  template<typename DataType>
  class CAMB {
  protected:
    using CoordinateType=tools::datatypes::strip_complex<DataType>;
    std::vector<CoordinateType> kcamb; //!< Wavenumbers read from CAMB file
    std::vector<std::vector<CoordinateType>> T; //!< Vector to store transfer functions:
    size_t cambCols[2] = {1,2}; //!< Columns of CAMB that we request:
    std::vector<tools::numerics::Interpolator<CoordinateType>> interpolator; //!< Interpolation functions:
    CoordinateType amplitude; //!< Amplitude of the initial power spectrum
    CoordinateType ns;        //!< tensor to scalar ratio of the initial power spectrum
    mutable CoordinateType kcamb_max_in_file; //!< Maximum CAMB wavenumber. If too small compared to grid resolution, Meszaros solution will be computed

  public:
    /*! \brief Vector which maps inputs to the correct transfer function.

       Elements are all zero if dmOnly = true,
       and just ascending if dmOnly = false. This prevents calls to operator() etc from going out of bounds if
       we switch off the use of extra transfer functions, and avoids having to tell all the fields to request
       a different transfer function in this case.*/
    std::vector<size_t> transferSwitch;
    const size_t nTransfers = 2; //!< Total number of transfer functions stored:
    bool dmOnly; //!< True if only using dark matter transfer function in the simulation
    bool dataRead; //!< True if we have read in data.

    //! Constructor. Assumes that only dark matter is being used, unless we specify otherwise.
    CAMB(bool useOnlyDM = true)
    {
        this->dmOnly = useOnlyDM;
        this->dataRead = false;
        // Initialise the interpolators and transfer functions vectors to be the
        // appropriate size, containing default constructed objects.
        if(useOnlyDM)
        {
            this->interpolator.resize(1);
            this->T.resize(1);
            transferSwitch.assign(this->nTransfers,0);
        }
        else
        {
            this->interpolator.resize(this->nTransfers);
            this->T.resize(this->nTransfers);
            for(size_t i = 0;i < this->T.size();i++)
            {
                transferSwitch.push_back(i);
            }
        }
    }

    //! Check that the supplied transfer function is valid (ie, that we don't request a transfer function which hasn't been imported)
    void ensureValidTransferID(const size_t transferID) const {
        if(transferSwitch[transferID] > this->T.size() - 1)
      {
        throw(std::runtime_error("Invalid transfer function requested."));
      }
    }

    //! Enables the use of the baryon transfer function
    void enableAllTransfers()
    {
        this->dmOnly = false;
        if(this->T.size() < 2)
        {
            this->kcamb.clear();
            this->T.clear();
        }
        this->T.resize(this->nTransfers);
        this->interpolator.resize(this->nTransfers);
        for(size_t i = 0;i < this->T.size();i++)
        {
            transferSwitch[i] = i;
        }
    }

    //! Checks whether the CAMB object actually stores any data:
    bool isUsable() const {
      return (kcamb.size() > 0);
    }

    //! Import data from CAMB file and initialise the interpolation functions used to compute the transfer functions:
    void read(std::string incamb, const CosmologicalParameters <CoordinateType> &cosmology) {
      readLinesFromCambOutput(incamb);
      for(size_t i = 0;i < this->T.size();i++)
      {
            this->interpolator[i].initialise(kcamb,this->T[i]);
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
      for(size_t i  =0;i < this->T.size();i++)
      {
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

      if(numCols == c_old_camb)
      {
        std::cerr << "Using pre 2015 CAMB transfer function" << std::endl;
      }
      else if(numCols == c_new_camb)
      {
        std::cerr << "Using post 2015 CAMB transfer function" << std::endl;
      }
      else
      {
        throw std::runtime_error("CAMB transfer file doesn't have a sensible number of rows and columns");
      }


      CoordinateType transferNormalisation = input[1]; // to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural
      // Copy file into vectors. Normalise both so that the DM tranfer function starts at 1.
      for (j = 0; j < input.size() / numCols; j++) {
        if (input[numCols * j] > 0) {
          // hard-coded to first two columns of CAMB file -
          kcamb.push_back(CoordinateType(input[numCols * j]));
          for(size_t i = 0;i < this->T.size();i++)
          {
            T[i].push_back(CoordinateType(input[numCols * j + this->cambCols[i]]) / transferNormalisation);
          }
        } else continue;
      }

      // Extend high-k range using Meszaros solution
      // This is a very naive approximation and a big warning will be issued if the power is actually evaluated
      // at these high k's (see operator() below).

      std::vector<CoordinateType> T_f;
      for(size_t i = 0;i < this->T.size();i++)
      {
        T_f.push_back(T[i].back());
      }
      kcamb_max_in_file = kcamb.back();
      CoordinateType keq = 0.01;
      while (kcamb.back() < 1e7) {
        kcamb.push_back(kcamb.back() + 1.0);
        CoordinateType kratio = kcamb.back() / kcamb_max_in_file;
        for(size_t i = 0;i < this->T.size();i++)
        {
            T[i].push_back(T_f[i] * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
        }
      }


    }

  public:

    //! \brief Evaluate power spectrum nTransfer at wavenumber k (Mpc/h), including the normalisation
    //! nTransfer = 0 -> DM
    //! nTransfer = 1 -> baryons.
    //! NB - for efficiency, we don't check that nTransfer is within bounds. This is handled by
    //! getPowerSpectrumForGrid, which does a single check before applying a guaranteed safe
    //! value to all fourier cells.
    CoordinateType operator()(CoordinateType k,int nTransfer = 0) const {
      CoordinateType linearTransfer;
      if (k != 0)
        linearTransfer = interpolator[transferSwitch[nTransfer]](k);
        // transferSwitch ensures that if extra transfer functions are switched off, then we always
        // call the dark matter interpolator, rather than falling out of bounds.
      else
        linearTransfer = 0.0;

      if (k > kcamb_max_in_file) {
        kcamb_max_in_file = std::numeric_limits<CoordinateType>().max();
      }

      return amplitude * powf(k, ns) * linearTransfer * linearTransfer;
    }

    //! Calculate and associate the theoretical power spectrum to a given grid
    std::shared_ptr<fields::Field<DataType, CoordinateType>>
    getPowerSpectrumForGrid(const grids::Grid<CoordinateType> &grid,size_t nTransfer = 0) const {
      /* Get the variance for each Fourier cell of the specified grid  */
      assert(kcamb.size() == T[0].size());
      // Check that the requested transfer function is actually available, and if it isn't, use the DM transfer function instead.
      // Do this here instead of inside operator() for efficiency.
      ensureValidTransferID(nTransfer);
      // transferSwitch ensures that if extra transfer functions are switched off, then we always
      // call the dark matter interpolator, rather than falling out of bounds.

      CoordinateType norm = getPowerSpectrumNormalizationForGrid(grid);

      auto P = std::make_shared<fields::Field<DataType, CoordinateType>>(grid, true);

      P->forEachFourierCell([norm, this,nTransfer]
                                (std::complex<CoordinateType>, CoordinateType kx, CoordinateType ky,
                                 CoordinateType kz) {
        CoordinateType k = sqrt(kx * kx + ky * ky + kz * kz);
        auto spec = std::complex<CoordinateType>((*this)(k,this->transferSwitch[nTransfer]) * norm, 0);
        return spec;
      });

      if (kcamb_max_in_file == std::numeric_limits<CoordinateType>().max()) {
        std::cerr << "WARNING: maximum k in CAMB input file is insufficient" << std::endl;
        std::cerr << "         extrapolating using naive Meszaros solution" << std::endl;
      }

      return P;

    }


  protected:
    //! Return the cosmology-dependent part of the normalisation of the power spectrum.
    void calculateOverallNormalization(const CosmologicalParameters <CoordinateType> &cosmology) {
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
    CoordinateType calculateLinearVarianceInSphere(CoordinateType radius,size_t nTrans = 0) const {
      // Bounds check
      ensureValidTransferID(nTrans);

      CoordinateType s = 0., k, t;

      CoordinateType amp = 9. / 2. / M_PI / M_PI;
      CoordinateType kmax = std::min(kcamb.back(), 200.0 / radius);
      CoordinateType kmin = kcamb[0];

      CoordinateType dk = (kmax - kmin) / 50000.;
      for (k = kmin; k < kmax; k += dk) {


      // Call the appropriate interpolation function for the requested transfer function
      // transferSwitch ensures that if extra transfer functions are switched off, then we always
      // call the dark matter interpolator, rather than falling out of bounds.
        t = interpolator[transferSwitch[nTrans]](k);

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
