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
    std::vector<CoordinateType> kcamb; /*!< Wavenumbers read from CAMB file */
    //std::vector<CoordinateType> Tcamb;  /*!< Transfer function at each k read from CAMB file, for DM field */
    //std::vector<CoordinateType> Tbar;/*!Transfer function for the baryons.*/
    //Vector to store transfer functions:
    std::vector<std::vector<CoordinateType>> T;
    //Columns of CAMB that we request:
    size_t cambCols[2] = {1,2};

    //Vector which maps inputs to the correct transfer function. Elements are all zero if dmOnly = true,
    //and just ascending if dmOnly = false. This prevents calls to operator() etc from going out of bounds if
    //we switch off the use of extra transfer functions, and avoids having to tell all the fields to request
    //a different transfer function in this case.
    std::vector<size_t> transferSwitch;
    //Interpolation functions:
    std::vector<toold::numerics::Interpolator<CoordinateType>> interpolator;
    //tools::numerics::Interpolator<CoordinateType> interpolator;
    //tools::numerics::Interpolator<CoordinateType> interpolator_baryons;
    CoordinateType amplitude; /*!< Amplitude of the initial power spectrum */
    CoordinateType ns;        /*!< tensor to scalar ratio of the initial power spectrum*/
    mutable CoordinateType kcamb_max_in_file; /*!< Maximum CAMB wavenumber. If too small compared to grid resolution, Meszaros solution will be computed */

  public:
    //Total number of transfer functions stored:
    const size_t nTransfers = 2;
    //Flags whether we want to only use the DM transfer function (true, default option) or to include baryons/others (false);
    bool dmOnly;

    //Constructor:
    CAMB(bool useOnlyDM = true)
    {
        this->dmOnly = useOnlyDM;
        //Initialise the interpolators and transfer functions vectors to be the
        //appropriate size, containing default constructed objects.
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

    //Function to switch to using both DM and baryons:
    void setDMOnly()
    {
        this->dmOnly = true;
        if(this->T.size() > 1)
        {
            std::cerr << "WARNING: imported CAMB data for additional transfer functions will be lost!!" << std::endl;
        }
        this->interpolator.resize(1);
        this->T.resize(1);
        transferSwitch.assign(this->T.size(),0);
    }
    //Function to enable using the baryon transfer function:
    void enableAllTransfers()
    {
        this->dmOnly = false;
        if(this->T.size() < 2)
        {
            std::cerr << "All transfer functions enabled. Remember to read in new CAMB data." << std::endl;
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

    //Checks whether the CAMB object actually stores any data:
    bool isUsable() const {
      return (kcamb.size() > 0);
    }

    //Import data from CAMB file and initialise the interpolation functions used to compute
    //the transfer functions:
    void read(std::string incamb, const CosmologicalParameters <CoordinateType> &cosmology) {
      readLinesFromCambOutput(incamb);
      for(int i = 0;i < this->T.size();i++)
      {
            this->interpolator[i].initialise(kcamb,this->T[i]);
      }
      //interpolator.initialise(kcamb, Tcamb);
      //interpolator_baryons.initialise(kcamb, Tbar);

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
  //This function imports data from a CAMB file, supplied as a file-name string argument (incamb). Both pre-2015 and post-2015
  //formats can be used, and the function will detect which.
    void readLinesFromCambOutput(std::string incamb) {
      kcamb.clear();
      for(int i  =0;i < this->T.size();i++)
      {
        T[i].clear();
      }


      // Dealing with the update of the CAMB TF output. Both are kept for backward compatibility.
      const int c_old_camb = 7; // number of columns in camb transfer function pre 2015
      const int c_new_camb = 13; // number of columns in camb transfer function post 2015
      int c;
      size_t j;

      //Import data from CAMB file:
      std::vector<double> input;

      io::getBuffer(input, incamb);

      //Check whether the input file is in the pre-2015 or post-2015 format (and throw an error if it is neither).
      if(input.size() > c_old_camb && input.size() % c_old_camb == 0){
        std::cerr << "Using pre 2015 CAMB transfer function" << std::endl;
        c = c_old_camb;
      } else if (input.size() > c_new_camb && input.size() % c_new_camb == 0){
        std::cerr << "Using post 2015 CAMB transfer function" << std::endl;
        c = c_new_camb;
      } else{
        throw std::runtime_error("CAMB transfer file doesn't have a sensible number of rows and columns");
      }


      CoordinateType ap = input[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural
      //Copy file into vectors. Normalise both so that the DM tranfer function starts at 1.
      for (j = 0; j < input.size() / c; j++) {
        if (input[c * j] > 0) {
          // hard-coded to first two columns of CAMB file -
          kcamb.push_back(CoordinateType(input[c * j]));
          for(int i = 0;i < this->T.size();i++)
          {
            T[i].push_back(CoordinateType(input[c * j + this->cambCols[i]]) / ap);
          }
          //Tcamb.push_back(CoordinateType(input[c * j + 1]) / ap);
          //Tbar.push_back(CoordinateType(input[c * j + 2]) / ap);
        } else continue;
      }

      // Extend high-k range using Meszaros solution
      // This is a very naive approximation and a big warning will be issued if the power is actually evaluated
      // at these high k's (see operator() below).

      std::vector<CoordinateType> T_f;
      for(int i = 0;i < this->T.size();i++)
      {
        T_f.push_back(T[i].back());
      }
      //CoordinateType Tcamb_f = Tcamb.back();
      //CoordinateType Tbar_f = Tbar.back();
      kcamb_max_in_file = kcamb.back();
      CoordinateType keq = 0.01;
      while (kcamb.back() < 1e7) {
        kcamb.push_back(kcamb.back() + 1.0);
        CoordinateType kratio = kcamb.back() / kcamb_max_in_file;
        for(int i = 0;i < this->T.size();i++)
        {
            T[i].push_back(T_f[i] * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
        }
        //Tcamb.push_back(Tcamb_f * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
        //Tbar.push_back(Tbar_f * pow(kratio, -2.0) * log(kcamb.back() / keq) / log(kcamb_max_in_file / keq));
      }


    }

  public:

    //! Evaluate power spectrum nTransfer at wavenumber k (Mpc/h), including the normalisation
    //! nTransfer = 0 -> DM
    //! nTransfer = 1 -> baryons.
    //!NB - for efficiency, we don't check that nTransfer is within bounds. This is handled by
    //! getPowerSpectrumForGrid, which does a single check before applying a guaranteed safe
    //! value to all fourier cells.
    CoordinateType operator()(CoordinateType k,int nTransfer = 0) const {
      CoordinateType linearTransfer;
      if (k != 0)
        linearTransfer = interpolator[transferSwitch[nTransfer]](k);
        //transferSwitch ensures that if extra transfer functions are switched off, then we always
        //call the dark matter interpolator, rather than falling out of bounds.
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
      //Check that the requested transfer function is actually available, and if it isn't, use the DM transfer function instead.
      //Do this here instead of inside operator() for efficiency.
      if(transferSwitch[nTransfer] > this->T.size() - 1)
      {
        std::cerr << "WARNING: invalid transfer function requested. Defaulting to dark matter transfer function." << std::endl;
        nTransfer = 0;
      }
      //transferSwitch ensures that if extra transfer functions are switched off, then we always
      //call the dark matter interpolator, rather than falling out of bounds.

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
    void calculateOverallNormalization(const CosmologicalParameters <CoordinateType> &cosmology) {
      CoordinateType ourGrowthFactor = growthFactor(cosmology);
      CoordinateType growthFactorNormalized = ourGrowthFactor / growthFactor(cosmologyAtRedshift(cosmology, 0));
      CoordinateType sigma8PreNormalization = calculateLinearVarianceInSphere(8.);
      CoordinateType linearRenormFactor = (cosmology.sigma8 / sigma8PreNormalization) * growthFactorNormalized;

      amplitude = linearRenormFactor * linearRenormFactor;

    }

  public:

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


    //Compute the power spectrum averaged over a spherical top hat window.
    CoordinateType calculateLinearVarianceInSphere(CoordinateType radius,size_t nTransfer = 0) const {
      //Bounds check
      if(transferSwitch[nTransfer] > this->T.size() - 1)
      {
        std::cerr << "WARNING: invalid transfer function requested. Defaulting to dark matter transfer function." << std::endl;
        nTransfer = 0;
      }

      CoordinateType s = 0., k, t;

      CoordinateType amp = 9. / 2. / M_PI / M_PI;
      CoordinateType kmax = std::min(kcamb.back(), 200.0 / radius);
      CoordinateType kmin = kcamb[0];

      CoordinateType dk = (kmax - kmin) / 50000.;
      for (k = kmin; k < kmax; k += dk) {


      //Call the appropriate interpolation function for the requested transfer function
      //transferSwitch ensures that if extra transfer functions are switched off, then we always
      //call the dark matter interpolator, rather than falling out of bounds.
        t = interpolator[transferSwitch[nTransfer]](k);

        //Multiply power spectrum by the fourier transform of the spherical top hat, to give the fourier transform
        //of the averaged (convolved) power spectrum over the sphere.
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
