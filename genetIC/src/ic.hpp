#ifndef _IC_HPP_INCLUDED
#define _IC_HPP_INCLUDED

#include <string>
#include <tuple>
#include <cassert>
#include <functional>
#include <algorithm>
#include <memory>
#include <limits>
#include <iostream>
#include <list>

#include "io/numpy.hpp"

#include <src/simulation/modifications/modificationmanager.hpp>

#include "src/tools/filesystem.h"

#include "src/simulation/field/randomfieldgenerator.hpp"
#include "src/simulation/field/multilevelfield.hpp"
#include "src/simulation/field/evaluator.hpp"

#include "src/simulation/particles/multilevelgenerator.hpp"
#include "src/simulation/particles/mapper/onelevelmapper.hpp"
#include "src/simulation/particles/mapper/twolevelmapper.hpp"
#include "src/simulation/particles/mapper/gasmapper.hpp"
#include "src/simulation/particles/mapper/graficmapper.hpp"
#include "src/simulation/particles/velocityoffsetgenerator.hpp"

#include "src/cosmology/camb.hpp"
#include "src/simulation/window.hpp"

using namespace std;

template<typename T>
class DummyICGenerator;


/*!
   \class ICGenerator
   \brief top level object responsible for coordinating the generation of initial conditions, including genetic modifications.

   This class exposes all methods accessible at user level through main.o

*/
template<typename GridDataType>
class ICGenerator {
protected:

  using T = tools::datatypes::strip_complex<GridDataType>;
  using GridPtrType = std::shared_ptr<grids::Grid<T>>;


  friend class DummyICGenerator<GridDataType>;

  cosmology::CosmologicalParameters<T> cosmology;
  multilevelcontext::MultiLevelContextInformation<GridDataType> multiLevelContext;

  //! Vector of output fields (defaults to just a single field)
  std::vector<std::shared_ptr<fields::OutputField<GridDataType>>> outputFields;
  int nFields;
  modifications::ModificationManager<GridDataType> modificationManager;
  std::vector<fields::RandomFieldGenerator<GridDataType>> randomFieldGenerator;

  //! Switch to re-direct calls to the baryon field to the DM field if we are not using it:
  std::vector<size_t> transferSwitch;



  cosmology::CAMB<GridDataType> spectrum;

  //! Velocity offset to be added uniformly to output (default 0,0,0)
  Coordinate<GridDataType> velOffset;

  //! Gadget particle types to be generated on each level (default 1)
  std::vector<unsigned int> gadgetTypesForLevels;

  //! If true, at output time assign a different gadget particle type to the flagged particles on the finest known level
  bool flaggedParticlesHaveDifferentGadgetType;
  unsigned int flaggedGadgetParticleType;

  //! DM supersampling to perform on zoom grid
  int supersample;

  //! Subsampling on base grid
  int subsample;

  //! Normalisation used for cell softening scale:
  T epsNorm = 0.01075; // Default value arbitrary to coincide with normal UW resolution


  io::OutputFormat outputFormat;
  string outputFolder, outputFilename;

  //! Track whether the random realisation has yet been made
  bool haveInitialisedRandomComponent;

  //! Enforce the exact power spectrum, as in Angulo & Pontzen 2016
  bool exactPowerSpectrum;

  /*!
   Stray" particles are high-res particles outside a high-res grid,
   constructed through interpolation of the surrounding low-res grid. Disabled by default.
   */
  bool allowStrayParticles;

  //! If true, the box are recentered on the last centered point in the parameter file
  // TODO This is currently a grafic only option and ought to be generic
  bool centerOnTargetRegion = false ;
  bool haveCentred = false;

  //! If true, then the code will generate baryons on all levels, rather than just the deepest level.
  bool baryonsOnAllLevels;

  //! Flags whether an output format has been specified or not, so that an error can be generated if not.
  bool formatSpecified = false;


  std::vector<size_t> flaggedParticles;
  std::vector<std::vector<size_t>> zoomParticleArray;


  //! Coordinates of the cell of current interest
  T x0 = 0.0, y0 = 0.0, z0 = 0.0; // Have to initialise these to something, otherwise it could cause unpredictable errors at runtime.

  //! Value of passive variable for refinement masks if needed
  T pvarValue = 1.0;

  //! High-pass filtering scale defined for variance calculations
  T variance_filterscale = -1.0;

  //! Numerical parameters for the algorithm of the quadratic modifications
  // These two parameters are for the quadratic algorithm (Rey and Pontzen 2017):
  // initial_number_steps encodes the number of steps always performed by the algorithm (similar to a burn-in)
  // 10 is a good trade-off to reach acceptable precision while avoiding extra steps.
  // precision is the target precision at which the quadratic modification will be achieved. More precision requires
  // more steps. In practice, we are limited by other errors in the code (namely recombination of low and high ks for fields)
  // of order of 1%. precision=0.1% is therefore sufficient and carries an acceptable numerical cost.
  int initial_number_steps = 10;
  T precision = 0.001;

  shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pMapper;
  shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pInputMapper;
  shared_ptr<multilevelcontext::MultiLevelContextInformation<GridDataType>> pInputMultiLevelContext;

  std::vector<shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>>> pParticleGenerator;

  using FieldType = fields::Field<GridDataType>;

  tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter;


public:
  //! \brief Main constructor for this class.
  //! Requires a single interpreter to initialise, from which it will receive input commands.
  ICGenerator(tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter) :
      modificationManager(multiLevelContext, cosmology, outputFields),
      pMapper(new particle::mapper::ParticleMapper<GridDataType>()),
      interpreter(interpreter) {
    // Set default parameters, in case input file didn't specify these:
    // By default, we assume there is only one field - the DM field:
    nFields = 1;
    outputFields.push_back(std::make_shared<fields::OutputField<GridDataType>>(multiLevelContext,0));
    // Link random field generator the dark matter field:
    randomFieldGenerator.emplace_back(*(outputFields[0]));
    transferSwitch.push_back(0); // Enables use of dark matter transfer function
    velOffset = {0,0,0};

    // By default, no input mapper is used:
    pInputMapper = nullptr;
    pInputMultiLevelContext = nullptr;

    // Default cosmological paramters:
    cosmology.hubble = 0.701;   // old default
    cosmology.OmegaBaryons0 = 0.0;
    cosmology.ns = 0.96;      // old default
    cosmology.TCMB = 2.725;
    haveInitialisedRandomComponent = false;

    // Default computational options:
    supersample = 1;
    subsample = 1;
    exactPowerSpectrum = false;
    allowStrayParticles = false;
    auto _pParticleGenerator = std::make_shared<particle::NullMultiLevelParticleGenerator<GridDataType>>();
    pParticleGenerator.push_back(_pParticleGenerator);
    flaggedParticlesHaveDifferentGadgetType = false;
    flaggedGadgetParticleType = 1;
    this->baryonsOnAllLevels = false; // Baryons only output on the finest level by default.
  }

  ~ICGenerator() {
  }

  //!\brief Sets the matter density fraction to in.
  void setOmegaM0(T in) {
    cosmology.OmegaM0 = in;
  }

  //!\brief Sets the temperature of the CMB today to in.
  void setTCMB(T in) {
    cosmology.TCMB = in;
  }

  //!\brief Sets the baryon density fraction to in.
  void setOmegaB0(T in) {
    cosmology.OmegaBaryons0 = in;
    // Now that we have gas, mapper may have changed:
    updateParticleMapper();
  }

  //!\brief Sets the dark energy density fraction to in.
  void setOmegaLambda0(T in) {
    cosmology.OmegaLambda0 = in;
  }

  //!\brief sets the present-day Hubble rate
  //! \param in = Hubble rate dimensionless parameter, h, where H0 = 100h kms^{-1}Mpc^{-1}.
  void setHubble(T in) {
    cosmology.hubble = in;
  }

  //! Enables the use of stray particles
  void setStraysOn() {
    allowStrayParticles = true;
  }

//!\brief Enables the use of a baryon density field.
  //! If flagged on, the code will extract the transfer function for baryons from CAMB separately from
  //! the dark matter transfer. This causes twice as many fields to be stored and generated, but
  //! produces more accurate results for baryons than assuming they follow the same transfer function
  //! (which holds only at late times).
  void setUsingBaryons() {
    if(this->nFields > 1) {
        std::cerr<< "Already using baryons!" << std::endl;
    }
    else {
        if(this->spectrum.dataRead) {
            throw(std::runtime_error("Cannot switch on baryons after transfer function data already read."));
        }
        this->nFields = 2;
        // Add the baryon field:
        outputFields.push_back(std::make_shared<fields::OutputField<GridDataType>>(multiLevelContext,1));
        randomFieldGenerator.emplace_back(*(outputFields[1]));
        transferSwitch.push_back(1);
        // Add a particle generator for the baryon field:
        auto _pParticleGenerator = std::make_shared<particle::NullMultiLevelParticleGenerator<GridDataType>>();
        pParticleGenerator.push_back(_pParticleGenerator);
        this->spectrum.enableAllTransfers();
    }
  }

  //! Enables outputting baryons on all levels, rather than only the deepest level.
  void setBaryonsOnAllLevels() {
    this->baryonsOnAllLevels = true;
  }
  //! Sets a whole-box velocity offset in km/s -- e.g. for testing AMR sensitivity to movement relative to grid structure
  void setVelocityOffset(T vx, T vy, T vz) {
    // convert to km a^1/2 s^-1 (gadget units internally in code)
    T conversion = pow(cosmology.scalefactor, -0.5);
    velOffset = {vx*conversion,vy*conversion,vz*conversion};
  }

  //! Sets the sigma8 parameter.
  void setSigma8(T in) {
    cosmology.sigma8 = in;
  }

  //! Set the normalisation of the cell softening scale. Defaults to 0.01075
  void setEpsNorm(T in) {
    this->epsNorm = in;
  }

  //! Add a higher resolution grid to the stack by supersampling the finest grid.
  /*! The power spectrum will not be taken into account in this grid
   * \param factor Factor by which the resolution must be increased compared to the finest grid
   */
  void setSupersample(int factor) {
    if (factor <= 0 ){
      throw std::runtime_error("Supersampling factor must be greater then zero");
    }
    supersample = factor;
    updateParticleMapper();
  }

  //! Add a lower resolution grid to the stack by subsampling the coarsest grid.
  /*! The power spectrum will not be taken into account in this grid
   * \param factor Factor by which the resolution will be downgraded compared to the coarsest grid
   */
  void setSubsample(int factor) {
    if (factor <= 0 ){
      throw std::runtime_error("Subsampling factor must be greater then zero");
    }

    subsample = factor;
    updateParticleMapper();
  }

  //! Defines the redshift at which initial conditions are generated.
  void setZ0(T z) {
    cosmology.redshift = z;
    cosmology.scalefactor = 1. / (cosmology.redshift + 1.);
  }

  //! If generating an extra IC file of a passive variable, sets its value.
  void setpvarValue(T value) {
    this->pvarValue = value;
  }

  //! Centres the grid on the currently selected region - used for grafic output.
  void setCenteringOnRegion(){
    this->centerOnTargetRegion = true;
    updateParticleMapper();
  }

  //! Sets the cutoff scale used in the variance filter
  void setVarianceFilteringScale(T filterscale){
    this->variance_filterscale = filterscale;
  }

  //! Define the base (coarsest) grid
  /*!
   * \param boxSize Physical size of box in Mpc
   * \param n Number of cells in the grid
   */
  virtual void initBaseGrid(T boxSize, size_t n) {
    assert(boxSize > 0);

    if (multiLevelContext.getNumLevels() > 0)
      throw std::runtime_error("Cannot re-initialize the base grid");


    if(haveInitialisedRandomComponent) {
        throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));
    }

    addLevelToContext(spectrum, boxSize, n);
    updateParticleMapper();

  }

  //! Sets ns, the scalar spectral index.
  void setns(T in) {
    cosmology.ns = in;
  }

  /*! Define the zoomed grid encompassing all flagged cells/particles
   * The zoom region is defined such that the flagged cells are roughly in the middle of it. Its physical size
   * must ensure that all flagged cells fit inside it.
   * \param zoomfac Ratio between the physical sizes of the base and zoom grid
   * \param n Number of cells in the zoom grid
   */
  void initZoomGrid(size_t zoomfac, size_t n) {
    if(haveInitialisedRandomComponent) {
        throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));
    }


    if (multiLevelContext.getNumLevels() < 1)
      throw std::runtime_error("Cannot initialise a zoom grid before initialising the base grid");

    grids::Grid<T> &gridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);
    int nAbove = (int) gridAbove.size;

    storeCurrentCellFlagsAsZoomMask(multiLevelContext.getNumLevels());
    vector<size_t> &newLevelZoomParticleArray = zoomParticleArray.back();

    if (newLevelZoomParticleArray.size() == 0)
      throw std::runtime_error("Cannot initialise zoom without marking particles to be replaced");

    // find boundaries
    Window<int> zoomWindow(gridAbove.getEffectiveSimulationSize(),
                           gridAbove.getCoordinateFromIndex(newLevelZoomParticleArray[0]));

    for (auto cell_id : newLevelZoomParticleArray) {
      zoomWindow.expandToInclude(gridAbove.getCoordinateFromIndex(cell_id));
    }

    size_t n_required = zoomWindow.getMaximumDimension();

    // Now see if the zoom the user chose is OK
    size_t n_user = nAbove / zoomfac;
    if (n_required > n_user && !allowStrayParticles) {
      throw (std::runtime_error(
          "Zoom particles do not fit in specified sub-box. Decrease zoom, or choose different particles"));
    }

    if (n_required < n_user)
      zoomWindow.expandSymmetricallyToSize(n_user);

    // The edges of zooms regions carry numerical errors due to interpolation between levels (see ref)
    // Do not use them if you can
    int borderSafety = 3;
    for (auto cell_id : zoomParticleArray.back()){
      if( ! zoomWindow.containsWithBorderSafety(gridAbove.getCoordinateFromIndex(cell_id), borderSafety)){
        std::cerr << "WARNING: Opening a zoom where flagged particles are within " << borderSafety <<
            " pixels of the edge. This is prone to numerical errors." << std::endl;
        break;
      }
    }



    auto lci = zoomWindow.getLowerCornerInclusive();
    initZoomGridWithLowLeftCornerAt(lci.x, lci.y, lci.z, zoomfac, n);

  }

  //! Gets the currently flagged cells on the specified level and converts them to a zoom mask.
  void storeCurrentCellFlagsAsZoomMask(size_t level) {
    assert(level > 0);

    if (zoomParticleArray.size() < level)
      zoomParticleArray.emplace_back();

    assert(zoomParticleArray.size() >= level);

    grids::Grid<T> &gridAbove = multiLevelContext.getGridForLevel(level - 1);

    vector<size_t> &levelZoomParticleArray = zoomParticleArray[level - 1];
    levelZoomParticleArray.clear();
    gridAbove.getFlaggedCells(levelZoomParticleArray);
  }

  /*! Define a zoomed grid with user defined coordinates
   * \param x0, y0, z0  Coordinates in pixel number of the lower left corner
   */
  void initZoomGridWithLowLeftCornerAt(int x0, int y0, int z0, size_t zoomfac, size_t n) {
    grids::Grid<T> &gridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);
    int nAbove = int(gridAbove.size);

    storeCurrentCellFlagsAsZoomMask(multiLevelContext.getNumLevels());

    vector<size_t> trimmedParticleArray;
    vector<size_t> &untrimmedParticleArray = zoomParticleArray.back();
    vector<size_t> zoomedParticleArray;

    int nCoarseCellsOfZoomGrid = nAbove / int(zoomfac);

    Coordinate<int> lowerCorner(x0, y0, z0);
    Coordinate<int> upperCornerExclusive = lowerCorner + nCoarseCellsOfZoomGrid;

    if (gridAbove.coversFullSimulation())
      upperCornerExclusive = gridAbove.wrapCoordinate(upperCornerExclusive);
    else {
      if (upperCornerExclusive.x > nAbove || upperCornerExclusive.y > nAbove || upperCornerExclusive.z > nAbove)
        throw std::runtime_error("Attempting to initialise a zoom grid that falls outside the parent grid");
    }


    size_t missed_particle = 0;

    Window<int> zoomWindow = Window<int>(gridAbove.getEffectiveSimulationSize(),
                                         lowerCorner, upperCornerExclusive);


    // Make list of the particles, excluding those that fall outside the new high-res box. Alternatively,
    // if allowStrayParticles is true, keep even those outside the high-res box but report the number
    // in this category.
    for (size_t i = 0; i < untrimmedParticleArray.size(); i++) {
      bool include = true;
      auto coord = gridAbove.getCoordinateFromIndex(untrimmedParticleArray[i]);
      if (!zoomWindow.contains(coord)) {
        missed_particle += 1;
        include = false;
      }

      if (include || allowStrayParticles) {
        trimmedParticleArray.push_back(untrimmedParticleArray[i]);
      }
    }

    if (missed_particle > 0) {
      cerr << "WARNING: the requested zoom particles do not all fit in the requested zoom window" << endl;
      if(allowStrayParticles) {
        cerr << "         of " << untrimmedParticleArray.size() << " particles, " << missed_particle << " will be interpolated from LR grid (stray particle mode)" << endl;
      } else {
        cerr << "         of " << untrimmedParticleArray.size() << " particles, " << missed_particle << " have been omitted" << endl;
      }

      cerr << "         to make a new zoom flag list of " << trimmedParticleArray.size() << endl;
    }

    zoomParticleArray.pop_back();
    zoomParticleArray.emplace_back(std::move(trimmedParticleArray));

    Coordinate<T> newOffsetLower = gridAbove.offsetLower + Coordinate<T>(x0, y0, z0) * gridAbove.cellSize;

    this->addLevelToContext(spectrum, gridAbove.thisGridSize / zoomfac, n, newOffsetLower);


    grids::Grid<T> &newGrid = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);

    cout << "Initialized a zoom region:" << endl;
    cout << "  Subbox length         = " << newGrid.thisGridSize << " Mpc/h" << endl;
    cout << "  n                     = " << newGrid.size << endl;
    cout << "  dx                    = " << newGrid.cellSize << endl;
    cout << "  Zoom factor           = " << zoomfac << endl;
    cout << "  Num particles         = " << untrimmedParticleArray.size() << endl;
    cout << "  Low-left corner in parent grid = " << lowerCorner << endl;
    cout << "  Low-left corner (h**-1 Mpc)    = " << newGrid.offsetLower.x << ", " << newGrid.offsetLower.y << ", "
         << newGrid.offsetLower.z << endl;

    updateParticleMapper();

    cout << "  Total particles = " << pMapper->size() << endl;

    // Update the cell flags. The goal is to flag
    // Flag all cells on the fine grid which are included in the zoom
    vector<size_t> transferredCells;
    gridAbove.makeProxyGridToMatch(newGrid)->getFlaggedCells(transferredCells);
    newGrid.flagCells(transferredCells);

  }

  //!\brief Adds a level to the multi-level context.
  /*!
  * \param spectrum - CAMB object used to handle import of transfer functions.
  * \param size - size of level in Mpc
  * \param nside - number of cells on one side of the level
  * \param offset - co-ordinate for the lower left hand corner of the grid.
  */
  virtual void
  addLevelToContext(const cosmology::CAMB<GridDataType> &spectrum, T size, size_t nside,
                    const Coordinate<T> &offset = {0, 0, 0}) {
    // This forwards to multiLevelContext but is required because it is overriden in DummyICGenerator,
    // which needs to ensure that grids are synchronised between two different contexts
    multiLevelContext.addLevel(spectrum, size, nside, offset);
    this->gadgetTypesForLevels.push_back(1);
  }

  //!\brief Set all fields to be randomised with the same seed.
  /*!
  \param - seed to be applied.
  */
  void setSeed(int seed) {
    for(size_t i = 0; i < outputFields.size(); i++) {
        this->setSeed(seed,i);
    }
  }
  //! Seed all fields in Fourier space with the same seed.
  void setSeedFourier(int seed) {
    for(size_t i = 0; i < outputFields.size(); i++) {
        this->setSeedFourier(seed,i);
    }
  }

  //! \brief Set the seed in real space.
  /*!
  \param seed - seed to set.
  \param nField - field to set the seed for. 0 = dark matter, 1 = baryons.
  */
  void setSeed(int seed,size_t nField) {
    if(nField > outputFields.size() - 1) {
        throw(std::runtime_error("Attempted to apply operations to field that has not been setup. Use baryon_tf_on to enable baryons."));
    }
    randomFieldGenerator[nField].seed(seed);
  }

  //! \brief Set the seed in fourier space
  /*!
  \param seed - seed to set.
  \param nField - field to set the seed for. 0 = dark matter, 1 = baryons.
  */
  void setSeedFourier(int seed,size_t nField) {
    if(nField > outputFields.size() - 1) {
        throw(std::runtime_error("Attempted to apply operations to field that has not been setup. Use baryon_tf_on to enable baryons."));
    }
    randomFieldGenerator[nField].seed(seed);
    randomFieldGenerator[nField].setDrawInFourierSpace(true);
    randomFieldGenerator[nField].setReverseRandomDrawOrder(false);
    randomFieldGenerator[nField].setParallel(false);
    }

  //! Enable parallel random field generation, at the cost of not being backwards-compatible
  void setSeedFourierParallel(int seed,size_t nField) {
    randomFieldGenerator[nField].seed(seed);
    randomFieldGenerator[nField].setDrawInFourierSpace(true);
    randomFieldGenerator[nField].setParallel(true);
    randomFieldGenerator[nField].setReverseRandomDrawOrder(false);
  }
  //! Apply only to the dark matter field:
  void setSeedFourierParallel(int seed) {
    for(size_t i = 0; i < this->outputFields.size(); i++) {
        this->setSeedFourierParallel(seed,i);
    }
  }

  //! Set the gadget particle type to be produced by the deepest level currently in the grid hiearchy
  void setGadgetParticleType(unsigned int type) {
    if(type>=6)
      throw(std::runtime_error("Gadget particle type must be between 0 and 5"));
    if(multiLevelContext.getNumLevels()==0)
      throw(std::runtime_error("Can't set a gadget particle type until a level has been initialised to set it for"));
    this->gadgetTypesForLevels[multiLevelContext.getNumLevels()-1] = type;
    this->updateParticleMapper();
  }

  //! Request that particles on the finest level are additionally split in gadget output such that
  //! flagged particles have a different particle type
  void setFlaggedGadgetParticleType(unsigned int type) {
    if(type>=6)
      throw(std::runtime_error("Gadget particle type must be between 0 and 5"));
    flaggedParticlesHaveDifferentGadgetType = true;
    flaggedGadgetParticleType = type;
    this->updateParticleMapper();
  }

//!\brief Specifies the seed in fourier space, but with reversed order of draws between real and imaginary part of complex numbers.
  /*!
   * Provided for compatibility problems as different compilers handle the draw order differently
   * \param seed - seed to set.
   * \param nField - field to set the seed for: 0 = Dark-Matter, 1 = Baryons.
   */
  void setSeedFourierReverseOrder(int seed,size_t nField) {
    if(nField > outputFields.size() - 1) {
        throw(std::runtime_error("Attempted to apply operations to field that has not been setup. Use baryon_tf_on to enable baryons."));
    }
    std::cerr << "*** Warning: seed_field_fourier_reverse and seedfourier_reverse are deprecated commands and should be avoided.***" << std::endl;
    randomFieldGenerator[nField].seed(seed);
    randomFieldGenerator[nField].setDrawInFourierSpace(true);
    randomFieldGenerator[nField].setReverseRandomDrawOrder(true);
    randomFieldGenerator[nField].setParallel(false);
  }

  //! Seed in reverse order for all fields with the same seed.
  void setSeedFourierReverseOrder(int seed) {
    for(size_t i = 0; i < outputFields.size(); i++) {
        this->setSeedFourierReverseOrder(seed,i);
    }
  }
  void setExactPowerSpectrumEnforcement() {
    exactPowerSpectrum = true;
  }

  //! \brief Obtain power spectrum from a CAMB data file
  /*!
  * \param cambFieldPath - string of path to CAMB file
  */
  void setCambDat(std::string cambFilePath) {
    spectrum.read(cambFilePath, cosmology);
  }

  //! Set the ouput directory to the supplied string
  void setOutDir(std::string outputPath) {
    outputFolder = outputPath;
  }

  //! Sets prefix for the name of all output files.
  void setOutName(std::string outputFilename_) {
    outputFilename = outputFilename_;
  }

  //! Set formats from handled formats in io namespace
  /*!
   * \param format  2 =  Gadget2, 3 = Gadget3, 4 = tipsy, 5 = grafic
   */
  void setOutputFormat(int format) {
    outputFormat = static_cast<io::OutputFormat>(format);
    this->formatSpecified = true;
    updateParticleMapper();
  }

  //! Returns the directory currently used for output.
  string getOutputPath() {
    ostringstream fname_stream;
    if (outputFilename.size() == 0) {
      fname_stream << outputFolder << "/IC_" << tools::datatypes::floatinfo<T>::name << "_z" << cosmology.redshift
                   << "_" << multiLevelContext.getGridForLevel(0).size;
    } else {
      fname_stream << outputFolder << "/" << outputFilename;
    }
    return fname_stream.str();
  }

  //! Throw an error if the level specified does not exist.
  void checkLevelExists(size_t level,size_t nField) {
    if(level > outputFields[nField]->getNumLevels() - 1) {
        throw(std::runtime_error("Specified level does not exist."));
    }
  }

  //! Throw an error if the code attempts to apply operations to the baryon field without having created it.
  void checkFieldExists(size_t nField) {
    if(nField > outputFields.size() - 1) {
        throw(std::runtime_error("Attempted to apply operations to field that has not been setup. Use baryon_tf_on to enable baryons."));
    }
  }

  //! Check if we have initialised the random component of all the fields - if we haven't, then initialise them.
  void initialiseRandomComponentIfUninitialised() {
    if (!haveInitialisedRandomComponent) {
        initialiseAllRandomComponents();
    }
  }

  //! Zeroes field values at a given level. Meant for debugging.
  virtual void zeroLevel(size_t level,size_t nField) {
    checkFieldExists(nField);
    cerr << "*** Warning: your script calls zeroLevel(" << level << "). This is intended for testing purposes only!"
         << endl;

    checkLevelExists(level,nField);
    initialiseRandomComponentIfUninitialised();

    auto &fieldData = outputFields[nField]->getFieldForLevel(level).getDataVector();
    std::fill(fieldData.begin(), fieldData.end(), 0);
  }

  //! For backwards compatibility with scripts that don't use baryon transfer function
  virtual void zeroLevel(size_t level)
  {
    this->zeroLevel(level,0);
  }

  //! \brief Imports random field data for a level from a supplies file.
  /*!
  * \param level - level on which to import the data
  * \param filename - string giving the path to the file to be imported.
  * \param nField - optional (default is field 0). Specifies which field to import the data into. 0 = dark matter, 1 = baryons.
  */
  virtual void importLevel(size_t level, std::string filename,size_t nField = 0) {
    checkFieldExists(nField);
    initialiseRandomComponentIfUninitialised();
    cerr << "Importing random field on level " << level << " from " <<filename << endl;
    checkLevelExists(level,nField);


    auto & levelField = outputFields[nField]->getFieldForLevel(level);
    levelField.loadGridData(filename);
    levelField.setFourier(false);
    cerr << "... success!" << endl;
  }

  //! For backwards compatibility with scripts that don't use baryon transfer function. Imports data into dark matter field.
  virtual void importLevel(size_t level,std::string filename) {
    this->importLevel(level,filename,0);
  }

  //! \brief Apply the power spectrum to a specific field.
  //! \param nField - field to apply power spectrum to. 0 = dark matter, 1 = baryons.
  virtual void applyPowerSpec(size_t nField) {
    checkFieldExists(nField);

    if (this->exactPowerSpectrum) {
      outputFields[nField]->enforceExactPowerSpectrum();
    } else {
      outputFields[nField]->applyPowerSpectrum();
    }
  }

  //! Applies appropriate power spectrum to all fields.
  virtual void applyPowerSpec() {
    for(size_t i = 0 ;i < this->outputFields.size(); i++) {
        this->applyPowerSpec(i);
    }
  }

  //! \brief Outputs data about a specified field to a file.
  /*!
  * \param level - level of multi-level context on which field is found.
  * \param TField - field to dump.
  * \param prefix - prefix to define filename of output files.
  */
  template<typename TField>
  void dumpGridData(size_t level, const TField &data,size_t nField,std::string prefix = "grid") {

    initialiseRandomComponentIfUninitialised();

    auto levelGrid = data.getGrid();

    // Output grid data:
    ostringstream filename;
    filename << outputFolder << "/" << prefix << "-" << level;
    if(nField > 0) {
        filename << "-field-" << nField;
    }
    filename << ".npy";

    data.dumpGridData(filename.str());


    // Output meta-data:
    filename.str("");
    filename << outputFolder << "/" << prefix << "-" << level << ".txt";
    filename << outputFolder << "/" << prefix << "-" << level;
    if(nField > 0) {
        filename << "-field-" << nField;
    }
    if(data.isFourier()) {
        filename << "-fourier";
    }
    filename << ".txt";

    ofstream ifile;
    ifile.open(filename.str());
    cerr << "Writing to " << filename.str() << endl;

    ifile << levelGrid.offsetLower.x << " " << levelGrid.offsetLower.y << " "
          << levelGrid.offsetLower.z << " " << levelGrid.thisGridSize << endl;
    ifile << "The line above contains information about grid level " << level << endl;
    ifile << "It gives the x-offset, y-offset and z-offset of the low-left corner and also the box length" << endl;
    ifile.close();
  }

  //! Overload that the dark matter only (for backwards compatibility):
  template<typename TField>
  void dumpGridData(size_t level, const TField& data) {
    this->dumpGridData(level,data,0);
  }


  //! Dumps field in a tipsy format
  /*!
  \param fname - file-name to use for data.
  \param nField - field to dump. 0 = dark matter, 1 = baryons.
  */
  virtual void saveTipsyArray(string fname,size_t nField = 0) {
    checkFieldExists(nField);
    io::tipsy::saveFieldTipsyArray(fname, *pMapper, *pParticleGenerator[nField], *(outputFields[nField]) );
  }

  //! For backwards compatibility. Dumpys dark matter in tipsy format.
  virtual void saveTipsyArray(std::string fname)
  {
    this->saveTipsyArray(fname,0);
  }

  //! Dumps field at a given level in a file named grid-level
  /*!
  * \param level - level of multi-level context to dump
  * \param nField - field to dump. 0 = dark matter, 1 = baryons.
  */
  virtual void dumpGrid(size_t level,size_t nField) {
    checkFieldExists(nField);
    checkLevelExists(level,nField);
    outputFields[nField]->toReal();
    dumpGridData(level, outputFields[nField]->getFieldForLevel(level),nField);
    outputFields[nField]->toFourier();
  }

  //! For backwards compatibility. Dumpts baryons to field at requested level to file named grid-level.
  virtual void dumpGrid(size_t level) {
    this->dumpGrid(level,0);
  }

  //! For backwards compatibility. Dumpts baryons to field at requested level to file named grid-level.
  virtual void dumpGridFourier(size_t level = 0) {
    this->dumpGridFourier(level,0);
  }

  //! Output the grid in Fourier space.
  virtual void dumpGridFourier(size_t level = 0,size_t nField = 0) {
    checkFieldExists(nField);
    checkLevelExists(level,nField);

    fields::Field<complex<T>, T> fieldToWrite = tools::numerics::fourier::getComplexFourierField(
        outputFields[nField]->getFieldForLevel(level));
    dumpGridData(level, fieldToWrite,nField);
  }

  //! Dumps power spectrum generated from the field and the theory at a given level in a .ps file
  /*!
  \param level - level of multi-level context to dump
  \param nField - field to dump. 0 = dark matter, 1 = baryons.
  */
  virtual void dumpPS(size_t level,size_t nField) {
    checkFieldExists(nField);
    checkLevelExists(level,nField);

    auto &field = outputFields[nField]->getFieldForLevel(level);
    field.toFourier();
    std::string filename;
    if(nField == 1) {
        filename = (getOutputPath() + "_" + ((char) (level + '0')) + "_baryons.ps");
    }
    else {
        filename = (getOutputPath() + "_" + ((char) (level + '0')) + ".ps");
    }
    cosmology::dumpPowerSpectrum(field,
                                 multiLevelContext.getCovariance(level,nField),
                                 filename.c_str());
  }

  //! For backwards compatibility. Dumps dark matter power spectrum on the requested level.
  virtual void dumpPS(size_t level) {
    this->dumpPS(level,0);
  }
  //! For backwards compatibility. Dumps dark matter power spectrum on level 0.
  virtual void dumpPS() {
    this->dumpPS(0,0);
  }

  //! Dumps mask information to numpy grid files
  virtual void dumpMask() {
    cerr << "Dumping mask grids" << endl;
    // this is ugly but it makes sure I can dump virtual grids if there are any.
    multilevelcontext::MultiLevelContextInformation<GridDataType> newcontext;
    if(this->centerOnTargetRegion){
      this->multiLevelContext.copyContextWithCenteredIntermediate(newcontext, Coordinate<T>(x0,y0,z0), 2,
                                                                  subsample, supersample);}
    else{
      this->multiLevelContext.copyContextWithCenteredIntermediate(newcontext, this->getBoxCentre(), 2,
                                                                  subsample, supersample);
    }

    auto dumpingMask = multilevelcontext::GraficMask<GridDataType, T>(&newcontext, this->zoomParticleArray);
    dumpingMask.calculateMask();


    auto maskfield = dumpingMask.convertToField();
    for(size_t level=0; level<newcontext.getNumLevels(); level++){
      dumpGridData(level, maskfield->getFieldForLevel(level),0,std::string("mask"));
    }
  }


  //! Initialises the particle generator for field nField
  /*!
  * \param nField - field to initialise particle generator for. 0 = dark matter, 1 = baryons.
  */
  virtual void initialiseParticleGenerator(size_t nField = 0) {
    checkFieldExists(nField);

    // in principle this could now be easily extended to slot in higher order PT or other
    // methods of generating the particles from the fields

    using GridLevelGeneratorType = particle::ZeldovichParticleGenerator<GridDataType>;
    using OffsetGeneratorType = particle::VelocityOffsetMultiLevelParticleGenerator<GridDataType>;

    pParticleGenerator[nField] = std::make_shared<
        particle::MultiLevelParticleGenerator<GridDataType, GridLevelGeneratorType>>( *(outputFields[nField]) , cosmology, epsNorm);
    if(velOffset!=Coordinate<GridDataType>({0,0,0})) {
      cerr << "Adding a velocity offset to the output" << endl;
      pParticleGenerator[nField] = std::make_shared<OffsetGeneratorType>(pParticleGenerator[nField], velOffset);
    }

  }

  //! Runs commands of a given parameter file to set up the input mapper
  void setInputMapper(std::string fname) {
    if (multiLevelContext.getNumLevels() == 0)
      throw std::runtime_error("Mapper relative command cannot be used before a basegrid has been initialised.");


    DummyICGenerator<GridDataType> pseudoICs(this);
    auto dispatch = interpreter.specify_instance(pseudoICs);
    ifstream inf;
    inf.open(fname);


    if (!inf.is_open())
      throw std::runtime_error("Cannot open IC paramfile for relative_to command");
    cerr << "******** Running commands in " << fname << " to work out relationship ***********" << endl;

    tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(fname));

    dispatch.run_loop(inf);
    cerr << *(pseudoICs.pMapper) << endl;
    cerr << "******** Finished with " << fname << " ***********" << endl;
    pInputMapper = pseudoICs.pMapper;
    pInputMultiLevelContext = std::make_shared<multilevelcontext::MultiLevelContextInformation<GridDataType>>
        (pseudoICs.multiLevelContext);
  }

  /*! Get the grid on which the output is defined for a particular level.
   *
   * This may differ from the grid on which the fields are defined either because there is an offset or
   * there are differences in the resolution between the output and the literal fields.
   */
  std::shared_ptr<grids::Grid<T>> getOutputGrid(int level = 0) {
    auto gridForOutput = multiLevelContext.getGridForLevel(level).shared_from_this();

    if (allowStrayParticles && level > 0) {
      gridForOutput = std::make_shared<grids::ResolutionMatchingGrid<T>>(gridForOutput,
                                                                         getOutputGrid(level - 1));
    }
    return gridForOutput;
  }

  //!\brief Reconstructs the particle mapper after changes have been made.
  /*!
  * For example - if baryons have been added in, a gas mapper needs to be added, and if the output format has been set or
  * changed, this needs to be accounted for. It should be called whenever such a change has been made.
  */
  void updateParticleMapper() {
    // TODO: This routine contains too much format-dependent logic and should be refactored so that the knowledge
    // resides somewhere in the io namespace



    size_t nLevels = multiLevelContext.getNumLevels();

    if (nLevels == 0)
      return;

    if (outputFormat == io::OutputFormat::grafic) {
      // Grafic format just writes out the grids in turn. Grafic mapper only center when writing grids.
      // All internal calculations are done with center kept constant at boxsize/2.
      pMapper = std::make_shared<particle::mapper::GraficMapper<GridDataType>>(multiLevelContext, Coordinate<T>(x0,y0,z0),
                                                                               subsample, supersample);
      return;
    }

    // make a basic mapper for the coarsest grid
    auto baseLevelMapper = std::shared_ptr<particle::mapper::OneLevelParticleMapper<GridDataType>>(
      new particle::mapper::OneLevelParticleMapper<GridDataType>(
        getOutputGrid(0)
      ));

    baseLevelMapper->setGadgetParticleType(this->gadgetTypesForLevels[0]);

    pMapper = baseLevelMapper;


    if (nLevels >= 2) {

      for (size_t level = 1; level < nLevels; level++) {

        auto pFine = std::shared_ptr<particle::mapper::OneLevelParticleMapper<GridDataType>>(
            new particle::mapper::OneLevelParticleMapper<GridDataType>(getOutputGrid(level)));

        pFine->setGadgetParticleType(this->gadgetTypesForLevels[level]);

        pMapper = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
            new particle::mapper::TwoLevelParticleMapper<GridDataType>(pMapper, pFine, zoomParticleArray[level - 1]));
      }
    }

    if(flaggedParticlesHaveDifferentGadgetType) {
      auto pFlagged = std::shared_ptr<particle::mapper::OneLevelParticleMapper<GridDataType>>(
        new particle::mapper::OneLevelParticleMapper<GridDataType>(getOutputGrid(nLevels-1)));

      pFlagged->setGadgetParticleType(flaggedGadgetParticleType);

      vector<size_t> flaggedCells;
      pMapper->getFinestGrid()->getFlaggedCells(flaggedCells);
      pMapper = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
        new particle::mapper::TwoLevelParticleMapper<GridDataType>(pMapper, pFlagged, flaggedCells ));
    }

    if (cosmology.OmegaBaryons0 > 0) {

      // Add gas only to the deepest level. Pass the whole pGrid
      // vector if you want to add gas to every level.

      typedef std::pair<std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>,
                        std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>> gasMapperType;

      gasMapperType gasMapper;
      if(this->baryonsOnAllLevels) {
        // If we want to output baryons on all levels:
        gasMapper = pMapper->addGas(cosmology.OmegaBaryons0 / (cosmology.OmegaM0),
                                       multiLevelContext.getAllGrids());
      }
      else {
        // Default
        gasMapper = pMapper->addGas(cosmology.OmegaBaryons0 / cosmology.OmegaM0,
                                       {multiLevelContext.getGridForLevel(nLevels - 1).shared_from_this()});
      }


      bool gasFirst = outputFormat == io::OutputFormat::tipsy;

      // graft the gas particles onto the start of the map
      if (gasFirst)
        pMapper = std::make_shared<particle::mapper::AddGasMapper<GridDataType>>(
            gasMapper.first, gasMapper.second, true);
      else
        pMapper = std::make_shared<particle::mapper::AddGasMapper<GridDataType>>(
            gasMapper.second, gasMapper.first, false);

    }

    // potentially resample the lowest-level DM grid. Again, this is theoretically
    // more flexible if you pass in other grid pointers.
    if (supersample > 1)
      pMapper = pMapper->superOrSubSampleDM(supersample,
                                            {multiLevelContext.getGridForLevel(nLevels - 1).shared_from_this()}, true);

    if (subsample > 1)
      pMapper = pMapper->superOrSubSampleDM(subsample, {multiLevelContext.getGridForLevel(0).shared_from_this()},
                                            false);

  }

  //! \brief Clear out existing flags in the particle mapper and reflag using the selected particles
  /*!
    If an input mapper has been specified, then we map the flags specified relative to that setup
    into flags specified in the current setup. This allows for seamless transition, with flags defined
    in one simulation setup to be interpreted in another setup.
  */
  void reflag() {

    if (pInputMapper != nullptr) {
      pMapper->unflagAllParticles();
      pInputMapper->unflagAllParticles();
      pInputMapper->flagParticles(flaggedParticles);
      pInputMapper->extendParticleListToUnreferencedGrids(multiLevelContext);
      pMapper->extendParticleListToUnreferencedGrids(*pInputMultiLevelContext);
    } else {
      pMapper->unflagAllParticles();
      pMapper->flagParticles(flaggedParticles);
    }
  }


  //! Transforms the grid field in particles and outputs them in the predefined format
  virtual void write() {
    using namespace io;

    initialiseRandomComponentIfUninitialised();
    for(size_t i = 0; i < outputFields.size(); i++) {
        initialiseParticleGenerator(i);
    }


    cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
    cerr << (*pMapper);

    T boxlen = multiLevelContext.getGridForLevel(0).periodicDomainSize;

    switch (outputFormat) {
      case OutputFormat::gadget2:
      case OutputFormat::gadget3:
        gadget::save<float>(getOutputPath() + ".gadget", boxlen, *pMapper,
                     pParticleGenerator,
                     cosmology, static_cast<int>(outputFormat));
        break;
      case OutputFormat::tipsy:
        tipsy::save(getOutputPath() + ".tipsy", boxlen, pParticleGenerator,
                    pMapper, cosmology);
        break;
      case OutputFormat::grafic:
        if(this->centerOnTargetRegion){
          std::cerr << "Replacing coarse grids with centered grids on " << Coordinate<T>(x0,y0,z0) <<  std::endl;
          grafic::save(getOutputPath() + ".grafic",
                       pParticleGenerator, multiLevelContext, cosmology, pvarValue, Coordinate<T>(x0,y0,z0),
                       subsample, supersample, zoomParticleArray,outputFields);
        } else {
          grafic::save(getOutputPath() + ".grafic",
                       pParticleGenerator, multiLevelContext, cosmology, pvarValue, this->getBoxCentre(),
                       subsample, supersample, zoomParticleArray,outputFields);
        }

        break;
      default:
        throw std::runtime_error("Unknown output format");
    }

  }

  //! \brief Generates white noise for the specified field, and then applies the appropriate power spectrum.
  /*!
  * \param nField - field to apply to. 0 = dark matter, 1 = baryons. Defaults to 0.
  */
  void initialiseRandomComponent(size_t nField = 0) {
    checkFieldExists(nField);
    randomFieldGenerator[nField].draw();
    applyPowerSpec(nField);


  }

  //! Initialise random components for all the fields.
  void initialiseAllRandomComponents()
  {
    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to re-draw the random field after it was already initialised"));


    for(size_t i = 0; i < outputFields.size(); i++) {
        this->initialiseRandomComponent(i);
    }
    haveInitialisedRandomComponent = true;
  }


protected:

  //! Returns the level of the highest resolution grid with particles currently flagged.
  size_t deepestLevelWithParticlesSelected() {
    return multiLevelContext.deepestLevelwithFlaggedCells();
  }

  //! Returns the number of levels in the current multi-level context.
  size_t deepestLevel() {
    return multiLevelContext.getNumLevels();
  }

  T get_wrapped_delta(T x0, T x1) {
    return multiLevelContext.getGridForLevel(0).getWrappedOffset(x0, x1);
  }


  //! Outputs the centre of the currently flagged region
  void getCentre() {
    x0 = 0;
    y0 = 0;
    z0 = 0;

    size_t level = deepestLevelWithParticlesSelected();
    Coordinate<T> centre = this->multiLevelContext.getGridForLevel(level).getFlaggedCellsCentre();
    x0 = centre.x;
    y0 = centre.y;
    z0 = centre.z;
    cerr << "Centre of region is " << setprecision(12) << x0 << " " << y0 << " " << z0 << endl;
  }


  //! Loads flagged particles from a file, erases and duplicates, and then flags these particles.
  /*!
  * Note that this function does not erase existing flags - for that purpose, loadParticleIdFile should
  * be called instead.
  */
  void appendParticleIdFile(std::string filename) {

    cerr << "Loading " << filename << endl;

    io::getBuffer(flaggedParticles, filename);
    size_t size = flaggedParticles.size();
    tools::sortAndEraseDuplicate(flaggedParticles);
    if (flaggedParticles.size() < size)
      cerr << "  ... erased " << size - flaggedParticles.size() << " duplicate particles" << endl;
    cerr << "  -> total number of particles is " << flaggedParticles.size() << endl;

    reflag();
  }

  //! Clears the currently flagged particles and then loads new flags from a file.
  void loadParticleIdFile(std::string filename) {
    flaggedParticles.clear();
    appendParticleIdFile(filename);
  }


public:


//! Load from a file new flagged particles
  void loadID(string fname) {
    loadParticleIdFile(fname);
    getCentre();
  }

  //! Append from a file new flagged particles
  void appendID(string fname) {
    appendParticleIdFile(fname);
    getCentre();
  }

  //! Output to a file the currently flagged particles
  virtual void dumpID(string fname) {
    std::vector<size_t> results;
    cerr << "dumpID using current mapper:" << endl;
    cerr << (*pMapper);
    pMapper->getFlaggedParticles(results);
    io::dumpBuffer(results, fname);
  }

  //! Defines the currently interesting coordinates using a particle ID
  void centreParticle(long id) {
    std::tie(x0, y0, z0) = multiLevelContext.getGridForLevel(0).getCentroidFromIndex(id);
    this->haveCentred = true;
  }

  //! Flag the nearest cell to the coordinates currently pointed at
  void selectNearest() {
    auto &grid = multiLevelContext.getGridForLevel(deepestLevel() - 1);
    pMapper->unflagAllParticles();
    size_t id = grid.getIndexFromPoint(Coordinate<T>(x0, y0, z0));
    cerr << "selectNearest " << x0 << " " << y0 << " " << z0 << " " << id << " " << endl;
    grid.flagCells({id});
  }

  //! Selects particles to flag according to a specified function of their co-ordinate.
  void select(std::function<bool(T, T, T)> inclusionFunction) {
    T delta_x, delta_y, delta_z;
    T xp, yp, zp;

    flaggedParticles.clear();


    // unflag all grids first. This can't be in the loop below in case there are subtle
    // relationships between grids (in particular the ResolutionMatchingGrid which actually
    // points to two levels simultaneously).
    for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level) {
      getOutputGrid(level)->unflagAllCells();
    }

    for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level) {
      std::vector<size_t> particleArray;
      auto grid = getOutputGrid(level);
      size_t N = grid->size3;
      for (size_t i = 0; i < N; i++) {
        std::tie(xp, yp, zp) = grid->getCentroidFromIndex(i);
        delta_x = get_wrapped_delta(xp, x0);
        delta_y = get_wrapped_delta(yp, y0);
        delta_z = get_wrapped_delta(zp, z0);
        if (inclusionFunction(delta_x, delta_y, delta_z))
          particleArray.push_back(i);
      }
      grid->flagCells(particleArray);

    }
  }

  //! Flag all cells contained in the sphere centered at the coordinates currently pointed at
  /*!
   * \param radius in Mpc/h
   * */
  void selectSphere(float radius) {
    T r2 = radius * radius;
    select([r2](T delta_x, T delta_y, T delta_z) -> bool {
      T r2_i = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
      return r2_i < r2;
    });

  }

  //! Flag all cells contained in the cube centered at the coordinates currently pointed at
  /*!
   * \param side in Mpc/h
   * */
  void selectCube(float side) {
    T side_by_2 = side / 2;
    select([side_by_2](T delta_x, T delta_y, T delta_z) -> bool {
      return abs(delta_x) < side_by_2 && abs(delta_y) < side_by_2 && abs(delta_z) < side_by_2;
    });
  }

  //! Expand the current flagged region by the specified number of cells
  void expandFlaggedRegion(size_t nCells) {
    if(multiLevelContext.getNumLevels()<1) {
      throw std::runtime_error("No grid has yet been initialised");
    }
    std::cerr << "Expand flagged region by "<<nCells<< " cells" << std::endl;
    for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level) {
      grids::Grid<T> &thisGrid = multiLevelContext.getGridForLevel(level);
      size_t initial_length = thisGrid.numFlaggedCells();
      thisGrid.expandFlaggedRegion(nCells);
      std::cerr << "  - level " << level << " increased number of flagged cells by " << thisGrid.numFlaggedCells() - initial_length
                << " (now " << thisGrid.numFlaggedCells() <<")" << std::endl;
    }

  }

  //! On simulations with more than one zoom level, adapt the upper level zooms to fit snuggly around the lower levels
  /*! The actual position of the zoom box is never moved; this routine just re-selects the zoom cells
   *
   * \param nPadCells - the number of cells padding to incorporate
   */
  void adaptMask(size_t nPadCells) {
    if(multiLevelContext.getNumLevels()<3) {
      throw std::runtime_error("Adapting the mask requires there to be more than one zoom level");
    }
    assert(this->zoomParticleArray.size()==multiLevelContext.getNumLevels()-1);
    for(size_t level = multiLevelContext.getNumLevels()-2; level>0; --level) {
      auto &sourceArray = this->zoomParticleArray[level];
      grids::Grid<T> &thisLevelGrid = multiLevelContext.getGridForLevel(level);
      grids::Grid<T> &aboveLevelGrid = multiLevelContext.getGridForLevel(level-1);

      // Create a minimal mask on the level above, consisting of all parent cells of the zoom cells on this level
      thisLevelGrid.unflagAllCells();
      thisLevelGrid.flagCells(sourceArray);
      auto proxyGrid = thisLevelGrid.makeProxyGridToMatch(aboveLevelGrid);

      std::vector<size_t> zoomParticleArrayLevelAbove;
      proxyGrid->getFlaggedCells(zoomParticleArrayLevelAbove);

      // zoomParticleArrayLevelAbove now has the minimal mask.
      // Expand the minimal mask by nPadCells in each direction
      aboveLevelGrid.unflagAllCells();
      aboveLevelGrid.flagCells(zoomParticleArrayLevelAbove);
      aboveLevelGrid.expandFlaggedRegion(nPadCells);
      zoomParticleArrayLevelAbove.clear();
      aboveLevelGrid.getFlaggedCells(zoomParticleArrayLevelAbove);
      this->zoomParticleArray[level-1]=zoomParticleArrayLevelAbove;
    }
    updateParticleMapper();
  }

  //! Define cell currently pointed at by coordinates
  void setCentre(T xin, T yin, T zin) {
    x0 = xin;
    y0 = yin;
    z0 = zin;
    this->haveCentred = true;
  }

  //! Returns the coordinate of the centre of the lowest resolution grid
  Coordinate<T> getBoxCentre(){
    T boxsize = multiLevelContext.getGridForLevel(0).thisGridSize;
    return Coordinate<T>(boxsize/2,boxsize/2,boxsize/2);
  }

  //! Calculates and prints chi^2 for the specified field.
  /*!
  * \param nField - field to compute chi^2 for. 0 = dark matter, 1 = baryons.
  */
  virtual void getFieldChi2(size_t nField){
    checkFieldExists(nField);

    initialiseRandomComponentIfUninitialised();

    T val = this->outputFields[nField]->getChi2();
    std::cerr << "Calculated chi^2 = " <<  val <<std::endl;
  }

  //! Print chi2 for all fields:
  virtual void getFieldChi2() {
    for(size_t i = 0; i < this->outputFields.size(); i++) {
        if(i == 1) {
            std::cout << "Baryons: "; // To differentiate baryon output from dark matter.
        }
        this->getFieldChi2(i);
    }
  }

  //! Calculate physical quantities of the field
  /*!
   * @param filterscale Filtering scale in Mpc if the quantity needs it
   * @param nField - field to calculate for (0 = dark matter [default], 1 = baryons)
   */
  void calculate(string name,size_t nField = 0) {
    checkFieldExists(nField);
    initialiseRandomComponentIfUninitialised();

    GridDataType val = modificationManager.calculateCurrentValueByName(name, this->variance_filterscale);
    cout << name << ": calculated value = " << val << endl;
  }

  //! Calculate unmodified property for the dark matter field:
  void calculate(std::string name) {
    // Don't really need this for the baryon field, as we aren't modifying it:
    this->calculate(name,0);
  }

  //! Define a modification to be applied to the field
  /*!
   * \param name String with the name of the modification to be added.
   * @param type  Modification can be relative to existing value or absolute
   * @param target Absolute target or factor by which the existing will be multiplied
   */
  virtual void modify(string name, string type, float target) {
    initialiseRandomComponentIfUninitialised();

    modificationManager.addModificationToList(name, type, target,this->initial_number_steps, this->precision, this->variance_filterscale);
  }

  //! Empty modification list
  void clearModifications() {
    modificationManager.clearModifications();
  }

  //! Apply the algorithm to produce the modified field
  virtual void applyModifications() {
    modificationManager.applyModifications(); // This handles propagating the modifications to the other fields too.
  }


  //! Apply the modifications, calculate the corresponding delta chi^2 and recombine low and high-ks between grids to write a particle output
  virtual void done() {

    //Should exit gracefully if no output format specified (will be checked later,
    //but since format variable is undefined if not specified, it might accidentally
    //produce something unpredictable):
    if(!this->formatSpecified)
    {
        throw(std::runtime_error("No output format specified!"));
    }

    initialiseRandomComponentIfUninitialised();

    if(modificationManager.hasModifications())
      applyModifications();

    std::cout << "\nOperations applied to all fields. Writing to disk...";

    write();
  }

  //TODO This method does not really belong here but in MultiLevelField class
  void reverse(size_t nField = 0) {
    checkFieldExists(nField);
    outputFields[nField]->reverse();
  }

  //!For backwards compatibility
  void reverse() {
    for(size_t i = 0; i < outputFields.size(); i++){
        this->reverse(i);
    }
  }

  //TODO What is this method doing? Looks like inverted initial conditions properties
  //STEPHEN - renaming because this function actually reseeds high k, not small k
  //void reseedSmallK(T kmax, int seed,size_t nField = 0) {
  //!Reseeds the high k modes, leaving low-k modes unchanged.
  void reseedHighK(T kmax, int seed,size_t nField = 0) {

    checkFieldExists(nField);

    T k2max = kmax * kmax;


    // take a copy of all the fields
    std::vector<FieldType> fieldCopies;
    for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level) {
      auto &field = outputFields[nField]->getFieldForLevel(level);
      field.toFourier();
      fieldCopies.emplace_back(field);
    }

    // remake the fields with the new seed
    randomFieldGenerator[nField].seed(seed);
    initialiseRandomComponent(nField);

    // copy back the old field
    for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level) {
      FieldType &oldField = fieldCopies[level];
      auto &field = outputFields[nField]->getFieldForLevel(level);
      int k2max_i = tools::getRatioAndAssertInteger(k2max, field.getGrid().getFourierKmin());
      field.forEachFourierCell([k2max_i, &oldField](std::complex<T> val, int kx, int ky, int kz) {
        int k2_i = kx * kx + ky * ky + kz * kz;
        if (k2_i < k2max_i && k2_i != 0) {
          val = oldField.getFourierCoefficient(kx, ky, kz);//This line actually restores the original seeds for the small k modes, not the high k
        }
        return val;
      });
    }

  }

  //! Reseed high k for all fields
  void reseedHighK(T kmax,int seed)
  {
    for(size_t i = 0; i < outputFields.size(); i++){
        this->reseedHighK(kmax,seed,i);
    }
  }

  //! Reverses the sign of the low-k modes.
  void reverseSmallK(T kmax,size_t nField = 0) {

    checkFieldExists(nField);

    T k2max = kmax * kmax;

    for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level) {
      auto &field = outputFields[nField]->getFieldForLevel(level);
      field.toFourier();

      field.forEachFourierCell([k2max](std::complex<T> val, T kx, T ky, T kz) {
        T k2 = kx * kx + ky * ky + kz * kz;
        if (k2 < k2max && k2 != 0) {
          val = -val;
        }
        return val;
      });

      cerr << "reverseSmallK: k reversal at " << sqrt(k2max) << endl;
    }
  }

  //! Reverse small k for all fields
  void reverseSmallK(T kmax){
    for(size_t i = 0; i < outputFields.size(); i++){
        this->reverseSmallK(kmax,i);
    }
  }
};

#endif
