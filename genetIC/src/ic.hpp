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

#include "tools/numerics/vectormath.hpp"
#include "tools/numerics/fourier.hpp"
#include "tools/filesystem.h"
#include "io/numpy.hpp"
#include "cosmology/parameters.hpp"
#include "cosmology/camb.hpp"
#include "simulation/window.hpp"
#include "simulation/particles/species.hpp"
#include "simulation/multilevelgrid/multilevelgrid.hpp"
#include "simulation/field/randomfieldgenerator.hpp"
#include "simulation/field/multilevelfield.hpp"
#include "simulation/field/evaluator.hpp"
#include "simulation/modifications/modificationmanager.hpp"
#include "simulation/particles/multilevelgenerator.hpp"
#include "simulation/particles/mapper/onelevelmapper.hpp"
#include "simulation/particles/mapper/twolevelmapper.hpp"
#include "simulation/particles/mapper/gasmapper.hpp"
#include "simulation/particles/mapper/graficmapper.hpp"
#include "simulation/particles/offsetgenerator.hpp"
#include "tools/logging.hpp"

using namespace std;

template<typename T>
class DummyICGenerator;


/*!
   \class ICGenerator
   \brief Top level object responsible for coordinating the generation of initial conditions, including genetic modifications.

   This class exposes all methods accessible at user level through main.o

*/
template<typename GridDataType>
class ICGenerator {
protected:

  using T = tools::datatypes::strip_complex<GridDataType>;
  using GridPtrType = std::shared_ptr<grids::Grid<T>>;


  friend class DummyICGenerator<GridDataType>;

  cosmology::CosmologicalParameters<T> cosmology; //!< Cosmological parameters
  multilevelgrid::MultiLevelGrid<GridDataType> multiLevelContext; //!< Class to keep track of the multi-level context of the simulation

  //! Vector of output fields (defaults to just a single field)
  std::vector<std::shared_ptr<fields::OutputField<GridDataType>>> outputFields;
  bool useBaryonTransferFunction{false}; //!< True if gas particles should use a different transfer function
  modifications::ModificationManager<GridDataType> modificationManager; //!< Handles applying modificaitons to the various fields.
  std::unique_ptr<fields::RandomFieldGenerator<GridDataType>> randomFieldGenerator; //!< Generate white noise for the output fields
  std::unique_ptr<cosmology::PowerSpectrum<GridDataType>> spectrum; //!< Transfer function data

  //! Velocity offset to be added uniformly to output (default 0,0,0)
  Coordinate<GridDataType> velOffset;

  //! Gadget particle types to be generated on each level (default 1)
  std::vector<unsigned int> gadgetTypesForLevels;

  //! If true, at output time assign a different gadget particle type to the flagged particles on the finest known level
  bool flaggedParticlesHaveDifferentGadgetType;

  //! Gadget type of flagged particles.
  unsigned int flaggedGadgetParticleType;

  //! DM supersampling to perform on deepest zoom grid
  int supersample = 1;

  //! Gas supersampling to perform on deepest zoom grid
  int supersampleGas = 1;

  //! Subsampling on base grid
  int subsample = 1;

  //! Normalisation used for cell softening scale:
  T epsNorm = 0.01075; // Default value arbitrary to coincide with normal UW resolution


  io::OutputFormat outputFormat = io::OutputFormat::unknown; //!< Output format used by the code for particle data.
  string outputFolder; //!< Name of folder for output files.
  string outputFilename; //!< Name of files for output.

  //! Track whether the random realisation has yet been made
  bool haveInitialisedRandomComponent;

  //! Enforce the exact power spectrum, as in Angulo & Pontzen 2016
  bool exactPowerSpectrum;

  /*! \brief True if stray particles are allowed.

   Stray" particles are high-res particles outside a high-res grid,
   constructed through interpolation of the surrounding low-res grid. Disabled by default.
   */
  bool allowStrayParticles;

  /*! \brief If true, the box are recentered on the last centered point in the parameter file

   Note the centre can be set either explicitly ('center' command) or implicitly (e.g. by loading a set of particles
   with 'IDfile' command)
  */
  bool centerOnTargetRegion = false;

  //! If true, then the code will generate baryons on all levels, rather than just the deepest level.
  bool baryonsOnAllLevels;


  /*! \brief If >0, the code will generate 'padding' regions around deep zooms.
   *
   * E.g. if you zoom by a factor of 4, a thin layer of particles will be added around just outside the zoom region,
   * which are zoomed by only a factor of 2. This helps mitigate against boundary effects. The value of autopad
   * specifies the number of cells to pad out by (at the intermediate resolution).
   */
  size_t autopad = 0;

  //! Zoom masks for each level.
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
  int initial_number_steps = 10; //!< Number of steps always used by quadratic modifications.
  T precision = 0.001; //!< Target precision required of quadratic modifications.

  //! Mapper that keep track of particles in the mulit-level context.
  shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pMapper = nullptr;
  //! Input mapper, used to relate particles in a different simulation to particles in this one.
  shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pInputMapper;
  //! Multi-level context of the input mapper.
  shared_ptr<multilevelgrid::MultiLevelGrid<GridDataType>> pInputMultiLevelContext;

  //! Particle generators for each particle species.
  particle::SpeciesToGeneratorMap<GridDataType> pParticleGenerator;

  using FieldType = fields::Field<GridDataType>;

  tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter; //!< Parser for parameter files


public:
  //! \brief Main constructor for this class.
  //! Requires a single interpreter to initialise, from which it will receive input commands.
  ICGenerator(tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter) :
    modificationManager(multiLevelContext, cosmology, nullptr),
    interpreter(interpreter) {

    // By default, we assume there is only one field - at first this will contain white noise:
    outputFields.push_back(
      std::make_shared<fields::OutputField<GridDataType>>(multiLevelContext, particle::species::whitenoise));

    // Link modifications to the white noise field
    modificationManager.bindToField(outputFields[0]);

    // Link random field generator to the white noise field
    randomFieldGenerator = std::make_unique<fields::RandomFieldGenerator<GridDataType>>(*outputFields[0]);

    velOffset = {0, 0, 0};

    // By default, no input mapper is used:
    pInputMapper = nullptr;
    pInputMultiLevelContext = nullptr;

    // Default cosmological paramters:
    cosmology.hubble = 0.701;   // old default
    cosmology.OmegaBaryons0 = 0.0;
    cosmology.ns = 0.96;      // old default
    cosmology.TCMB = 2.725;
    cosmology.scalefactorAtDecoupling = 1./150.; // NB - this is NOT photon decoupling, it's baryon decoupling
    cosmology.OmegaM0 = 0.3;
    outputFolder = "./";

    // only an approximate value is needed to initialise a gas temperature
    
    haveInitialisedRandomComponent = false;

    // Default computational options:
    exactPowerSpectrum = false;
    allowStrayParticles = false;
    flaggedParticlesHaveDifferentGadgetType = false;
    flaggedGadgetParticleType = 1;
    this->baryonsOnAllLevels = false; // Baryons only output on the finest level by default.

  }

  //! Destructor
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
    multiLevelContext.allowStrays = true;
    allowStrayParticles = true;
  }

//!\brief Enables the use of a baryon density field.
  //! If flagged on, the code will extract the transfer function for baryons from CAMB separately from
  //! the dark matter transfer. This causes twice as many fields to be stored and generated, but
  //! produces more accurate results for baryons than assuming they follow the same transfer function
  //! (which holds only at late times).
  void setUsingBaryons() {
    this->useBaryonTransferFunction = true;
  }

  //! Enable autopad (i.e. generation of thin layers of intermediate resolution particles around zoom regions)
  void setAutopad(size_t nCells) {
    this->autopad = nCells;
    this->updateParticleMapper();
  }

  //! Enables outputting baryons on all levels, rather than only the deepest level.
  void setBaryonsOnAllLevels() {
    this->baryonsOnAllLevels = true;
  }

  //! Sets a whole-box velocity offset in km/s -- e.g. for testing AMR sensitivity to movement relative to grid structure
  void setOffset(T vx, T vy, T vz) {
    // convert to km a^1/2 s^-1 (gadget units internally in code)
    T conversion = pow(cosmology.scalefactor, -0.5);
    velOffset = {vx * conversion, vy * conversion, vz * conversion};
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
    if (factor <= 1) {
      throw std::runtime_error("Supersampling factor must be greater then one");
    }
    supersample = factor;
    updateParticleMapper();
  }

  //! Add a higher resolution grid to the stack by supersampling the finest grid.
  /*! The power spectrum will not be taken into account in this grid
   * \param factor Factor by which the resolution must be increased compared to the finest grid
   */
  void setSupersampleGas(int factor) {
    if (factor <= 1) {
      throw std::runtime_error("Supersampling factor must be greater then one");
    }
    supersampleGas = factor;
    updateParticleMapper();
  }

  //! Add a lower resolution grid to the stack by subsampling the coarsest grid.
  /*! The power spectrum will not be taken into account in this grid
   * \param factor Factor by which the resolution will be downgraded compared to the coarsest grid
   */
  void setSubsample(int factor) {
    if (factor <= 0) {
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
  void setCenteringOnRegion() {
    this->centerOnTargetRegion = true;
    updateParticleMapper();
  }

  //! Sets the cutoff scale used in the variance filter
  void setVarianceFilteringScale(T filterscale) {
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


    if (haveInitialisedRandomComponent) {
      throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));
    }

    logging::entry() << "Initialized the base grid:" << endl;
    logging::entry() << "  Box length            = " << boxSize << " Mpc/h" << endl;
    logging::entry() << "  n                     = " << n << endl;
    logging::entry() << "  dx                    = " << boxSize/n << endl;

    addLevelToContext(boxSize, n);
    updateParticleMapper();

  }

  //! Sets ns, the scalar spectral index.
  void setns(T in) {
    cosmology.ns = in;
  }

  /*! \brief Define the zoomed grid encompassing all flagged cells/particles.

   * The zoom region is defined such that the flagged cells are roughly in the middle of it. Its physical size
   * must ensure that all flagged cells fit inside it.
   * \param zoomfac Ratio between the physical sizes of the base and zoom grid
   * \param n Number of cells in the zoom grid
   */
  void initZoomGrid(size_t zoomfac, size_t n) {
    if (haveInitialisedRandomComponent) {
      throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));
    }


    if (multiLevelContext.getNumLevels() < 1)
      throw std::runtime_error("Cannot initialise a zoom grid before initialising the base grid");

    grids::Grid<T> &actualGridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);
    grids::Grid<T> &outputGridAbove = multiLevelContext.getOutputGridForLevel(multiLevelContext.getNumLevels() - 1);

    // Note everything that follows assumes the previous level's output grid (relative to which we must specify the
    // cells to zoom) has the same resolution as the actual grid (relative to which we specify the new grid).
    assert(actualGridAbove.cellSize == outputGridAbove.cellSize);

    int nAbove = (int) actualGridAbove.size;

    storeCurrentCellFlagsAsZoomMask(multiLevelContext.getNumLevels());
    vector<size_t> &newLevelZoomParticleArray = zoomParticleArray.back();

    if (newLevelZoomParticleArray.size() == 0)
      throw std::runtime_error("Cannot initialise zoom without marking particles to be replaced");

    // find boundaries
    Window<int> zoomWindow(outputGridAbove.getEffectiveSimulationSize(),
                           outputGridAbove.getCoordinateFromIndex(newLevelZoomParticleArray[0]));

    for (auto cell_id : newLevelZoomParticleArray) {
      zoomWindow.expandToInclude(outputGridAbove.getCoordinateFromIndex(cell_id));
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
    for (auto cell_id : zoomParticleArray.back()) {
      if (!zoomWindow.containsWithBorderSafety(outputGridAbove.getCoordinateFromIndex(cell_id), borderSafety)) {
        logging::entry(logging::level::warning) << "WARNING: Opening a zoom where flagged particles are within " << borderSafety <<
                  " pixels of the edge. This is prone to numerical errors." << std::endl;
        break;
      }
    }


    auto lci = zoomWindow.getLowerCornerInclusive();

    // The zoom grid is opened on top of the previous actual grid, not the output grid,
    // so we now need to convert. We set the domain size to be the grid above, so that we
    // can ensure the new zoom window is fully contained within that (unless the grid above
    // is the base level, in which case the wrapping is already done and no clamping should
    // be applied).

    Window<int> windowOnActualGrid(actualGridAbove.size,
      actualGridAbove.getCoordinateFromPoint(outputGridAbove.getCentroidFromCoordinate(lci)), lci+n_user);

    if(!outputGridAbove.coversFullSimulation())
      windowOnActualGrid.clampToFundamentalDomain(); // ensure we're not trying to wrap around when we can't!

    lci = windowOnActualGrid.getLowerCornerInclusive(); // fetch the updated lower corner coordinates

    initZoomGridWithLowLeftCornerAt(lci.x, lci.y, lci.z, zoomfac, n);

  }

  //! Gets the currently flagged cells on the specified level and converts them to a zoom mask.
  void storeCurrentCellFlagsAsZoomMask(size_t level) {
    assert(level > 0);

    if (zoomParticleArray.size() < level)
      zoomParticleArray.emplace_back();

    assert(zoomParticleArray.size() >= level);

    grids::Grid<T> &gridAbove = multiLevelContext.getOutputGridForLevel(level - 1);

    vector<size_t> &levelZoomParticleArray = zoomParticleArray[level - 1];
    levelZoomParticleArray.clear();
    gridAbove.getFlaggedCells(levelZoomParticleArray);
  }

  /*! \brief Define a zoomed grid with user defined coordinates
   * \param x0, y0, z0  Coordinates of the lower left corner in terms of the calculation grid from the level above
   * \param zoomfac The factor by which the zoom window is smaller than the level above
   * \param n The number of cells per side in the new zoom window
   */
  void initZoomGridWithLowLeftCornerAt(int x0, int y0, int z0, size_t zoomfac, size_t n) {
    grids::Grid<T> &calculationGridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);
    grids::Grid<T> &outputGridAbove = multiLevelContext.getOutputGridForLevel(multiLevelContext.getNumLevels() - 1);

    int nAbove = int(calculationGridAbove.size);

    storeCurrentCellFlagsAsZoomMask(multiLevelContext.getNumLevels());

    vector<size_t> trimmedParticleArray;
    vector<size_t> &untrimmedParticleArray = zoomParticleArray.back();
    vector<size_t> zoomedParticleArray;

    int nCoarseCellsOfZoomGrid = nAbove / int(zoomfac);

    Coordinate<int> lowerCorner(x0, y0, z0);
    Coordinate<int> upperCornerExclusive = lowerCorner + nCoarseCellsOfZoomGrid;

    if (x0 < 0 || y0 < 0 || z0 < 0)
      throw std::runtime_error("The newly initialised zoom level is not fully contained within the previous level");

    if (calculationGridAbove.coversFullSimulation())
      upperCornerExclusive = calculationGridAbove.wrapCoordinate(upperCornerExclusive-1)+1; // wrap the last included cell
    else {
      if (upperCornerExclusive.x > nAbove || upperCornerExclusive.y > nAbove || upperCornerExclusive.z > nAbove)
        throw std::runtime_error("The newly initialised zoom level is not fully contained within the previous level");
      // Note even if "strays" are enabled, this is considered illegal since it doesn't make sense to be putting high resolution
      // information on top of interpolated information from >1 level above.
    }

    size_t missed_particle = 0;

    T EPSILON = calculationGridAbove.cellSize*1e-6;

    Window<T> zoomWindow = Window<T>(calculationGridAbove.periodicDomainSize,
                                     calculationGridAbove.getCentroidFromCoordinate(lowerCorner) - EPSILON,
                                     calculationGridAbove.getCentroidFromCoordinate(upperCornerExclusive-1) + EPSILON);


    // Make list of the particles, excluding those that fall outside the new high-res box. Alternatively,
    // if allowStrayParticles is true, keep even those outside the high-res box but report the number
    // in this category.
    for (size_t i = 0; i < untrimmedParticleArray.size(); i++) {
      bool include = true;
      auto coord = outputGridAbove.getCentroidFromIndex(untrimmedParticleArray[i]);
      if (!zoomWindow.contains(coord)) {
        missed_particle += 1;
        include = false;
      }

      if (include || allowStrayParticles) {
        trimmedParticleArray.push_back(untrimmedParticleArray[i]);
      }
    }

    if (missed_particle > 0) {
      logging::entry(logging::level::warning) << "WARNING: the requested zoom particles do not all fit in the requested zoom window" << endl;
      if (allowStrayParticles) {
        logging::entry() << "         of " << untrimmedParticleArray.size() << " particles, " << missed_particle
             << " will be interpolated from LR grid (stray particle mode)" << endl;
      } else {
        logging::entry() << "         of " << untrimmedParticleArray.size() << " particles, " << missed_particle
             << " have been omitted" << endl;
      }

      logging::entry() << "         to make a new zoom flag list of " << trimmedParticleArray.size() << endl;
    }

    zoomParticleArray.pop_back();
    zoomParticleArray.emplace_back(std::move(trimmedParticleArray));

    Coordinate<T> newOffsetLower =
      calculationGridAbove.offsetLower + Coordinate<T>(x0, y0, z0) * calculationGridAbove.cellSize;

    this->addLevelToContext(calculationGridAbove.thisGridSize / zoomfac, n, newOffsetLower);


    grids::Grid<T> &newGrid = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);

    logging::entry() << "Initialized a zoom region:" << endl;
    logging::entry() << "  Subbox length         = " << newGrid.thisGridSize << " Mpc/h" << endl;
    logging::entry() << "  n                     = " << newGrid.size << endl;
    logging::entry() << "  dx                    = " << newGrid.cellSize << endl;
    logging::entry() << "  Zoom factor           = " << zoomfac << endl;
    logging::entry() << "  Num particles         = " << untrimmedParticleArray.size() << endl;
    logging::entry() << "  Low-left corner in parent grid = " << lowerCorner << endl;
    logging::entry() << "  Low-left corner (h**-1 Mpc)    = " << newGrid.offsetLower.x << ", " << newGrid.offsetLower.y << ", "
         << newGrid.offsetLower.z << endl;
    updateParticleMapper();

    logging::entry() << "  Total particles       = " << pMapper->size() << endl;

    // Flag all cells on the fine grid which were included in the zoom
    vector<size_t> transferredCells;
    calculationGridAbove.makeProxyGridToMatch(newGrid)->getFlaggedCells(transferredCells);
    newGrid.flagCells(transferredCells);

  }

  //!\brief Adds a level to the multi-level context.
  /*!
  * \param size - size of level in Mpc
  * \param nside - number of cells on one side of the level
  * \param offset - co-ordinate for the lower left hand corner of the grid.
  */
  virtual void
  addLevelToContext(T size, size_t nside, const Coordinate<T> &offset = {0, 0, 0}) {
    // This forwards to multiLevelContext but is required because it is overriden in DummyICGenerator,
    // which needs to ensure that grids are synchronised between two different contexts
    multiLevelContext.addLevel(size, nside, offset);
    this->gadgetTypesForLevels.push_back(1);
  }

  //! \brief Set up the random field generator to work in real space with a specified random seed
  /*!
  \param seed - seed to use
  */
  void setSeed(int seed) {
    randomFieldGenerator->seed(seed);
    randomFieldGenerator->setDrawInFourierSpace(false);
    randomFieldGenerator->setReverseRandomDrawOrder(false);
    randomFieldGenerator->setParallel(false);
  }


  //! \brief Set up the random field generator to work in Fourier space
  /*!
  \param seed - seed to use
  */
  void setSeedFourier(int seed) {
    randomFieldGenerator->seed(seed);
    randomFieldGenerator->setDrawInFourierSpace(true);
    randomFieldGenerator->setReverseRandomDrawOrder(false);
    randomFieldGenerator->setParallel(false);
  }


  //! \brief Set up the random field generator to work in Fourier space, using a parallel algorithm with a specified random seed
  /*!
  Note this is faster than using the older setSeedFourier setup, but not backwards- compatible
  (i.e. the resulting field will be different)
  \param seed - seed to use
  */
  void setSeedFourierParallel(int seed) {
    randomFieldGenerator->seed(seed);
    randomFieldGenerator->setDrawInFourierSpace(true);
    randomFieldGenerator->setParallel(true);
    randomFieldGenerator->setReverseRandomDrawOrder(false);
  }

  //!\brief Specifies the seed in fourier space, but with reversed order of draws between real and imaginary part of complex numbers.
  /*!
   * Provided for historical compatibility issues (where we discovered too late an initialization ambiguity in C++)
   * \param seed - seed to set.
   */
  void setSeedFourierReverseOrder(int seed) {
    std::cerr
      << "*** Warning: seed_field_fourier_reverse and seedfourier_reverse are deprecated commands and should be avoided.***"
      << std::endl;
    randomFieldGenerator->seed(seed);
    randomFieldGenerator->setDrawInFourierSpace(true);
    randomFieldGenerator->setReverseRandomDrawOrder(true);
    randomFieldGenerator->setParallel(false);
  }

  //! Enables exact power spectrum enforcement.
  void setExactPowerSpectrumEnforcement() {
    exactPowerSpectrum = true;
  }


  //! Set the gadget particle type to be produced by the deepest level currently in the grid hierarchy
  void setGadgetParticleType(unsigned int type) {
    if (type >= 6)
      throw (std::runtime_error("Gadget particle type must be between 0 and 5"));
    if (multiLevelContext.getNumLevels() == 0)
      throw (std::runtime_error("Can't set a gadget particle type until a level has been initialised to set it for"));
    this->gadgetTypesForLevels[multiLevelContext.getNumLevels() - 1] = type;
    this->updateParticleMapper();
  }

  //! Request that particles on the finest level are additionally split in gadget output such that
  //! flagged particles have a different particle type
  void setFlaggedGadgetParticleType(unsigned int type) {
    if (type >= 6)
      throw (std::runtime_error("Gadget particle type must be between 0 and 5"));
    flaggedParticlesHaveDifferentGadgetType = true;
    flaggedGadgetParticleType = type;
    this->updateParticleMapper();
  }






  //! \brief Obtain power spectrum from a CAMB data file
  /*!
  * \param cambFieldPath - string of path to CAMB file
  */
  void setCambDat(std::string cambFilePath) {
    spectrum = std::make_unique<cosmology::CAMB<GridDataType>>(this->cosmology, cambFilePath);
    this->multiLevelContext.setPowerspectrumGenerator(*spectrum);
  }

  //! \brief Use a power law spectrum with specified amplitude
  void setPowerLawAmplitude(T amplitude) {
    spectrum = std::make_unique<cosmology::PowerLawPowerSpectrum<GridDataType>>(this->cosmology, amplitude);
    this->multiLevelContext.setPowerspectrumGenerator(*spectrum);
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
  void setOutputFormat(io::OutputFormat format) {
    outputFormat = format;
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
  void checkLevelExists(size_t level, size_t nField) {
    if (level > outputFields[nField]->getNumLevels() - 1) {
      throw (std::runtime_error("Specified level does not exist."));
    }
  }

  //! Throw an error if the code attempts to apply operations to the baryon field without having created it.
  void checkFieldExists(size_t nField) {
    if (nField > outputFields.size() - 1) {
      throw (std::runtime_error(
        "Attempted to apply operations to field that has not been setup. Use baryon_tf_on to enable baryons."));
    }
  }

  //! Check if we have initialised the random component of all the fields - if we haven't, then initialise them.
  void initialiseRandomComponentIfUninitialised() {
    if (!haveInitialisedRandomComponent) {
      initialiseAllRandomComponents();
    }
  }

  //! Zeroes field values at a given level. Meant for debugging.
  virtual void zeroLevel(size_t level, size_t nField) {
    checkFieldExists(nField);
    logging::entry() << "*** Warning: your script calls zeroLevel(" << level << "). This is intended for testing purposes only!"
         << endl;

    checkLevelExists(level, nField);
    initialiseRandomComponentIfUninitialised();

    auto &fieldData = outputFields[nField]->getFieldForLevel(level).getDataVector();
    std::fill(fieldData.begin(), fieldData.end(), 0);
  }

  //! For backwards compatibility with scripts that don't use baryon transfer function
  virtual void zeroLevel(size_t level) {
    this->zeroLevel(level, 0);
  }

  //! \brief Imports overdensity field data for a given level from a supplied file.
  /*!
  * \param level - level on which to import the data
  * \param filename - string giving the path to the file to be imported.
  */
  virtual void importLevel(size_t level, std::string filename) {

    initialiseRandomComponentIfUninitialised();
    logging::entry() << "Importing random field on level " << level << " from " << filename << endl;
    checkLevelExists(level, 0);

    auto &levelField = outputFields[0]->getFieldForLevel(level);
    levelField.loadGridData(filename);
    levelField.setFourier(false);
    levelField.toFourier();
    outputFields[0]->applyTransferRatioOneLevel(particle::species::dm, particle::species::whitenoise, level); // transform back to whitenoise
    logging::entry() << "... success!" << endl;
  }


  //! Applies appropriate power spectrum to all fields.
  virtual void applyPowerSpec() {
    // This only makes sense if our output field so far is white noise
    assert(outputFields[0]->getTransferType() == particle::species::whitenoise);
    assert(outputFields.size()==1);

    if(useBaryonTransferFunction) {
      // make a copy of the white noise field:
      outputFields.emplace_back(std::make_shared<fields::OutputField<GridDataType>>(*outputFields[0]));
      outputFields[1]->applyPowerSpectrumFor(particle::species::baryon);
      outputFields[0]->applyPowerSpectrumFor(particle::species::dm);
    } else {
      outputFields[0]->applyPowerSpectrumFor(particle::species::all);
    }

    // Check everything has worked as expected (e.g. no accidental pointer copying)
    if(useBaryonTransferFunction)
      assert(outputFields[1]->getTransferType() == particle::species::baryon);

  }

  //! \brief Outputs data about a specified field to a file.
  /*!
  * \param level - level of multi-level context on which field is found.
  * \param TField - field to dump.
  * \param prefix - prefix to define filename of output files.
  */
  template<typename TField>
  void dumpGridData(size_t level, const TField &data, std::string prefix = "grid") {

    initialiseRandomComponentIfUninitialised();

    auto levelGrid = data.getGrid();

    // Output grid data:
    ostringstream filename;
    filename << outputFolder << "/" << prefix << "-" << level;
    filename << ".npy";

    data.dumpGridData(filename.str());


    // Output meta-data:
    filename.str("");
    // filename << outputFolder << "/" << prefix << "-" << level << ".txt";
    filename << outputFolder << "/" << prefix << "-info-" << level;
    if (data.isFourier()) {
      filename << "-fourier";
    }
    filename << ".txt";

    ofstream ifile;
    ifile.open(filename.str());

    ifile << levelGrid.offsetLower.x << " " << levelGrid.offsetLower.y << " "
          << levelGrid.offsetLower.z << " " << levelGrid.thisGridSize << endl;
    ifile << "The line above contains information about grid level " << level << endl;
    ifile << "It gives the x-offset, y-offset and z-offset of the low-left corner and also the box length" << endl;
    ifile.close();
  }


  //! Dumps overdensity of specified field in a tipsy format
  /*!
  \param fname - file-name to use for data.
  \param nField - field to dump. 0 = dark matter, 1 = baryons.
  */
  virtual void saveTipsyArray(string fname, size_t nField = 0) {
    checkFieldExists(nField);
    io::tipsy::saveFieldTipsyArray(fname, *pMapper, *pParticleGenerator[particle::dm], *(outputFields[nField]));
  }

  //! Dumps dark matter overdensity in tipsy format.
  virtual void saveTipsyArray(std::string fname) {
    this->saveTipsyArray(fname, 0);
  }

  //! Dumps overdensity field for a given species at a given level to file
  /*!
  * \param level - level in the grid hierarchy to dump
  * \param species - the field to dump
  */
  virtual void dumpGrid(size_t level, particle::species species ) {
    auto & field = this->getOutputFieldForSpecies(species);
    field.toReal();
    dumpGridData(level, field.getFieldForLevel(level));
  }

  //! Dumps dark matter overdensity field at a given level to file
  /*!
  * \param level - level in the grid hierarchy to dump
  */
  virtual void dumpGrid(size_t level) {
    dumpGrid(level, particle::species::dm);
  }

  virtual void dumpVelocityX(size_t level) {
    dumpGridData(level, *(this->pParticleGenerator[particle::species::dm]->getGeneratorForLevel(
      level).getGeneratedFields()[0]), "vx");
  }

  //! For backwards compatibility. Dumpts baryons to field at requested level to file named grid-level.
  virtual void dumpGridFourier(size_t level = 0) {
    this->dumpGridFourier(level, particle::species::dm);
  }

  //! Output the grid in Fourier space.
  virtual void dumpGridFourier(size_t level, particle::species species) {
    auto & field = this->getOutputFieldForSpecies(species);
    field.toFourier();
    fields::Field<complex<T>, T> fieldToWrite = tools::numerics::fourier::getComplexFourierField(field.getFieldForLevel(level));
    dumpGridData(level, fieldToWrite);
  }

  //! Dumps power spectrum generated from the field and the theory at a given level in a .ps file
  /*!
  \param level - level of multi-level context to dump
  \param species - the type of particle (in case of multiple transfer functions)
  */
  virtual void dumpPS(size_t level, particle::species species) {
    //int nField = static_cast<size_t>(species);
    //checkLevelExists(level, nField);

    auto &field = getOutputFieldForSpecies(species).getFieldForLevel(level);
    field.toFourier();
    std::string filename;
    if (species == particle::species::baryon) {
      filename = (getOutputPath() + "_" + ((char) (level + '0')) + "_baryons.ps");
    } else {
      filename = (getOutputPath() + "_" + ((char) (level + '0')) + ".ps");
    }

    if(!useBaryonTransferFunction)
      species = particle::species::all;

    cosmology::dumpPowerSpectrum(field,
                                 *multiLevelContext.getCovariance(level, species),
                                 filename.c_str());
  }

  //! For backwards compatibility. Dumps dark matter power spectrum on the requested level.
  virtual void dumpPS(size_t level) {
    this->dumpPS(level, particle::species::dm);
  }

  //! For backwards compatibility. Dumps dark matter power spectrum on level 0.
  virtual void dumpPS() {
    this->dumpPS(0, particle::species::dm);
  }

  //! Dumps mask information to numpy grid files
  virtual void dumpMask() {
    // this is ugly but it makes sure I can dump virtual grids if there are any.
    multilevelgrid::MultiLevelGrid<GridDataType> newcontext;
    if (this->centerOnTargetRegion) {
      this->multiLevelContext.copyContextWithCenteredIntermediate(newcontext, Coordinate<T>(x0, y0, z0), 2,
                                                                  subsample, supersample);
    } else {
      this->multiLevelContext.copyContextWithCenteredIntermediate(newcontext, this->getBoxCentre(), 2,
                                                                  subsample, supersample);
    }

    auto dumpingMask = multilevelgrid::GraficMask<GridDataType, T>(&newcontext, this->zoomParticleArray);
    dumpingMask.calculateMask();


    auto maskfield = dumpingMask.convertToField();
    for (size_t level = 0; level < newcontext.getNumLevels(); level++) {
      dumpGridData(level, maskfield->getFieldForLevel(level), "mask");
    }
  }


  //! Initialises the particle generator for a given species, connecting it to one of the output fields
  /*!
  * \param species - species of particle
  */
  virtual void initialiseParticleGenerator(particle::species species) {
    auto & outputField = getOutputFieldForSpecies(species);

    if(outputField.getTransferType()!=species) {
      // This species (probably baryons) actually just uses the output field for a different particle species
      // (probably DM).
      //
      // We don't want multiple particle generators based on the same output field -- that would be very wasteful!
      //      // So, just take a copy of the particle generator pointer.
      ensureParticleGeneratorInitialised(outputField.getTransferType());
      pParticleGenerator[species] = pParticleGenerator[outputField.getTransferType()];
    } else {

      initialiseParticleGeneratorUsingField(species, outputField);

    }
  }

  virtual void initialiseParticleGeneratorUsingField(particle::species species, fields::OutputField<GridDataType> & field) {
    // in principle this could now be easily extended to slot in higher order PT or other
    // methods of generating the particles from the fields

    using GridLevelGeneratorType = particle::ZeldovichParticleGenerator<GridDataType>;
    using OffsetGeneratorType = particle::OffsetMultiLevelParticleGenerator<GridDataType>;

    pParticleGenerator[species] = std::make_shared<
      particle::MultiLevelParticleGenerator<GridDataType, GridLevelGeneratorType>>(field, cosmology, epsNorm);

    Coordinate<GridDataType> posOffset;

    if (this->centerOnTargetRegion) {
      posOffset = {-x0, -y0, -z0};
      posOffset += getBoxCentre();
      posOffset = getOutputGrid(0)->wrapPoint(posOffset);
    }

    const Coordinate<GridDataType> zeroCoord;

    if (posOffset != zeroCoord || velOffset != zeroCoord) {
      pParticleGenerator[species] = std::make_shared<OffsetGeneratorType>(pParticleGenerator[species], posOffset,
                                                                          velOffset);
    }
  }

  virtual void ensureParticleGeneratorInitialised(particle::species species) {
    if (pParticleGenerator.count(species)==0) {
      initialiseParticleGenerator(species);
    }
  }

  fields::OutputField<GridDataType> & getOutputFieldForSpecies(particle::species species) {
    if(!useBaryonTransferFunction)
      species = particle::species::all;

    for(auto outputField : this->outputFields) {
      if(outputField->getTransferType()==species)
        return *outputField;
    }
    throw(std::out_of_range("Cannot find an output field for species"));
  }


  //! Initialises the particle generators for all species, if not already done
  virtual void ensureParticleGeneratorInitialised() {
    assert(outputFields.size()>0);
    ensureParticleGeneratorInitialised(particle::dm);

    // If there is an appropriate field available, we use it for baryon output. If not, the baryon particle generator
    // will just point back to the DM particle generator.
    ensureParticleGeneratorInitialised(particle::baryon);

    if(!useBaryonTransferFunction)
      assert(pParticleGenerator[particle::baryon]==pParticleGenerator[particle::dm]);

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
    logging::entry() << endl;
    logging::entry() << "+ Input mapper: computing geometry from " << fname << endl;
    {
      tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(fname));
      logging::IndentWhileInScope temporaryIndent;
      logging::entry() << endl;
      dispatch.run_loop(inf);
    }
    logging::entry() << endl;


#ifdef DEBUG_INFO
    logging::entry() << endl;
    logging::entry() << "Input mapper retrieved:" << endl;
    logging::entry() << *(pseudoICs.pMapper);
    logging::entry() << endl;
#endif


    pInputMapper = pseudoICs.pMapper;
    pInputMultiLevelContext = std::make_shared<multilevelgrid::MultiLevelGrid<GridDataType>>
      (pseudoICs.multiLevelContext);
  }

  /*! \brief Get the grid on which the output is defined for a particular level.
   *
   * This may differ from the grid on which the fields are defined either because there is an offset or
   * there are differences in the resolution between the output and the literal fields.
   */
  std::shared_ptr<grids::Grid<T>> getOutputGrid(int level = 0) {
    return multiLevelContext.getOutputGridForLevel(level).shared_from_this();
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
      pMapper = std::make_shared<particle::mapper::GraficMapper<GridDataType>>(multiLevelContext,
                                                                               Coordinate<T>(x0, y0, z0),
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

        auto pOutputGridThisLevel = getOutputGrid(level);

        auto pFine = std::shared_ptr<particle::mapper::OneLevelParticleMapper<GridDataType>>(
          new particle::mapper::OneLevelParticleMapper<GridDataType>(pOutputGridThisLevel));

        pFine->setGadgetParticleType(this->gadgetTypesForLevels[level]);

        pMapper = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
          new particle::mapper::TwoLevelParticleMapper<GridDataType>(pMapper, pFine, zoomParticleArray[level - 1]));
      }
    }

    if (flaggedParticlesHaveDifferentGadgetType) {
      auto pFlagged = std::shared_ptr<particle::mapper::OneLevelParticleMapper<GridDataType>>(
        new particle::mapper::OneLevelParticleMapper<GridDataType>(getOutputGrid(nLevels - 1)));

      pFlagged->setGadgetParticleType(flaggedGadgetParticleType);

      vector<size_t> flaggedCells;
      pMapper->getFinestGrid()->getFlaggedCells(flaggedCells);
      pMapper = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
        new particle::mapper::TwoLevelParticleMapper<GridDataType>(pMapper, pFlagged, flaggedCells));
    }

    decltype(pMapper) gasMapper, dmMapper;

    if (cosmology.OmegaBaryons0 > 0) {
      decltype(multiLevelContext.getAllGrids()) gridsToAddBaryonsTo;

      if (this->baryonsOnAllLevels)
        gridsToAddBaryonsTo = multiLevelContext.getAllGrids();
      else
        gridsToAddBaryonsTo = {multiLevelContext.getGridForLevel(nLevels - 1).shared_from_this()};
      // Default - only add gas on the finest level

      std::tie(gasMapper, dmMapper) = pMapper->splitMass(cosmology.OmegaBaryons0 / (cosmology.OmegaM0),
                                                         gridsToAddBaryonsTo);

      gasMapper = superOrSubSample(gasMapper, supersampleGas, 1);
      dmMapper = superOrSubSample(dmMapper, supersample, subsample);

      bool gasFirst = outputFormat == io::OutputFormat::tipsy;

      // graft the gas particles onto the start (or end) of the map
      if (gasFirst)
        pMapper = std::make_shared<particle::mapper::AddGasMapper<GridDataType>>(
          gasMapper, dmMapper, true);
      else
        pMapper = std::make_shared<particle::mapper::AddGasMapper<GridDataType>>(
          dmMapper, gasMapper, false);

    } else {
      pMapper = superOrSubSample(pMapper, supersample, subsample);
    }

    if (autopad > 0) {
      pMapper = pMapper->insertIntermediateResolutionPadding(2, autopad);
    }

  }

  auto superOrSubSample(decltype(pMapper) pForSubmap, int supersample, int subsample) {
    if (supersample > 1)
      pForSubmap = pForSubmap->superOrSubSample(supersample,
                                                {multiLevelContext.getGridForLevel(
                                                  multiLevelContext.getNumLevels() - 1).shared_from_this()}, true);

    if (subsample > 1)
      pForSubmap = pForSubmap->superOrSubSample(subsample, {multiLevelContext.getGridForLevel(0).shared_from_this()},
                                                false);

    return pForSubmap;
  }

  //! \brief Clear flags on all cells on all grids
  void clearCellFlags() {
    pMapper->unflagAllParticles();

    if (pInputMapper != nullptr)
      pInputMapper->unflagAllParticles();

  }

  //! \brief Flag cells which correspond to the particle IDs given
  /*!
    If an input mapper has been specified, then we map the flags specified relative to that setup
    into flags specified in the current setup. This allows for flags defined in one simulation setup
    to be interpreted in another setup.
  */
  void flagCellsCorrespondingToParticles(const std::vector<size_t> & flaggedParticles) {
    if (pInputMapper != nullptr) {
      pInputMapper->flagParticles(flaggedParticles);
      pInputMapper->extendParticleListToUnreferencedGrids(multiLevelContext);
      pMapper->extendParticleListToUnreferencedGrids(*pInputMultiLevelContext);
    } else {
      pMapper->flagParticles(flaggedParticles);
    }
  }


  //! Outputs the ICs in the defined format, creating appropriate particle generators if required
  virtual void write() {
    using namespace io;

    initialiseRandomComponentIfUninitialised();
    applyPowerSpec();
    ensureParticleGeneratorInitialised();

    logging::entry() << "Writing output; number dm particles=" << pMapper->size_dm()
         << ", number gas particles=" << pMapper->size_gas() << endl;
#ifdef DEBUG_INFO
    logging::entry() << (*pMapper);
#endif

    T boxlen = multiLevelContext.getGridForLevel(0).periodicDomainSize;
    Coordinate<T> centre;

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
        centre = this->getBoxCentre();

        if (this->centerOnTargetRegion) {
          // Unlike particle formats, for grafic we need to move the fields *on the grid* rather than just
          // at the point of particle generation.
          centre = Coordinate<T>(x0, y0, z0);
        }

        grafic::save(getOutputPath() + ".grafic",
                     pParticleGenerator, multiLevelContext, cosmology, pvarValue, centre,
                     subsample, supersample, zoomParticleArray, outputFields);
        break;
      default:
        throw std::runtime_error("Unknown output format");
    }

    logging::entry() << "Finished writing initial conditions" << endl;

  }

  //! Initialise random components for all the fields.
  virtual void initialiseAllRandomComponents() {
    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to re-draw the random field after it was already initialised"));


    // Draw the white noise field (writes to outputFields[0])
    randomFieldGenerator->draw();

    // Enforce exact unit variance if requested
    if(exactPowerSpectrum)
      outputFields[0]->enforceExactPowerSpectrum();

    // Make copies of the field ready for additional transfer functions (such as baryons)
    for (size_t i = 1; i < outputFields.size(); i++) {
      outputFields[i]->copyData(*(outputFields[0]));
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

  //! Gets the difference between the arguments, wrapped using the periodicity scale.
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
    logging::entry() << "Centre of region is " << setprecision(12) << x0 << " " << y0 << " " << z0 << endl;
  }


  //! Loads flagged particles from a file, erases and duplicates, and then flags these particles.
  /*!
  * Note that this function does not erase existing flags - for that purpose, loadParticleIdFile should
  * be called instead.
  */
  void appendParticleIdFile(std::string filename) {
#ifdef DEBUG_INFO
    logging::entry() << "Loading " << filename << endl;
#endif
    std::vector<size_t> flaggedParticles;
    io::getBuffer(flaggedParticles, filename);
    size_t size = flaggedParticles.size();
    tools::sortAndEraseDuplicate(flaggedParticles);

    flagCellsCorrespondingToParticles(flaggedParticles);
#ifdef DEBUG_INFO
    logging::entry() << "  -> total number of particle IDs loaded is " << flaggedParticles.size() << endl;
    flaggedParticles.clear();
    pMapper->getFlaggedParticles(flaggedParticles);
    logging::entry() << "     total number of accessible cells flagged is " << flaggedParticles.size() << endl;
#endif
  }

  //! Clears the currently flagged particles and then loads new flags from a file.
  void loadParticleIdFile(std::string filename) {
    clearCellFlags();
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
#ifdef DEBUG_INFO
    logging::entry() << "dumpID using current mapper" << endl;
    logging::entry() << (*pMapper);
#endif
    pMapper->getFlaggedParticles(results);
    io::dumpBuffer(results, fname);
  }

  //! Defines the currently interesting coordinates using a particle ID
  void centreParticle(long id) {
    std::tie(x0, y0, z0) = multiLevelContext.getGridForLevel(0).getCentroidFromIndex(id);
  }

  //! Flag the nearest cell to the coordinates currently pointed at
  void selectNearest() {
    auto &grid = multiLevelContext.getGridForLevel(deepestLevel() - 1);
    pMapper->unflagAllParticles();
    size_t id = grid.getIndexFromPoint(Coordinate<T>(x0, y0, z0));
#ifdef DEBUG_INFO
    logging::entry() << "selectNearest " << x0 << " " << y0 << " " << z0 << " " << id << " " << endl;
#endif
    grid.flagCells({id});
  }

  //! Selects particles to flag according to a specified function of their co-ordinate.
  void select(std::function<bool(T, T, T)> inclusionFunction) {
    T delta_x, delta_y, delta_z;
    T xp, yp, zp;

    // unflag all grids first. This can't be in the loop below in case there are subtle
    // relationships between grids (in particular the ResolutionMatchingGrid which actually
    // points to two levels simultaneously).
    for (size_t level = 0; level < multiLevelContext.getNumLevels(); ++level) {
      getOutputGrid(level)->unflagAllCells();
    }

    for (size_t level = 0; level < multiLevelContext.getNumLevels(); ++level) {
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
    if (multiLevelContext.getNumLevels() < 1) {
      throw std::runtime_error("No grid has yet been initialised");
    }
    logging::entry() << "Expand flagged region by " << nCells << " cells" << std::endl;
    for (size_t level = 0; level < multiLevelContext.getNumLevels(); ++level) {
      grids::Grid<T> &thisGrid = multiLevelContext.getGridForLevel(level);
      size_t initial_length = thisGrid.numFlaggedCells();
      thisGrid.expandFlaggedRegion(nCells);
      logging::entry() << "  - level " << level << " increased number of flagged cells by "
                << thisGrid.numFlaggedCells() - initial_length
                << " (now " << thisGrid.numFlaggedCells() << ")" << std::endl;
    }

  }

  //! On simulations with more than one zoom level, adapt the upper level zooms to fit snuggly around the lower levels
  /*! The actual position of the zoom box is never moved; this routine just re-selects the zoom cells
   *
   * \param nPadCells - the number of cells padding to incorporate
   */
  void adaptMask(size_t nPadCells) {
    if (multiLevelContext.getNumLevels() < 3) {
      throw std::runtime_error("Adapting the mask requires there to be more than one zoom level");
    }
    assert(this->zoomParticleArray.size() == multiLevelContext.getNumLevels() - 1);
    for (size_t level = multiLevelContext.getNumLevels() - 2; level > 0; --level) {
      auto &sourceArray = this->zoomParticleArray[level];
      auto thisLevelGrid = multiLevelContext.getOutputGridForLevel(level).withIndependentFlags();
      auto aboveLevelGrid = multiLevelContext.getOutputGridForLevel(level - 1).withIndependentFlags();

      // Create a minimal mask on the level above, consisting of all parent cells of the zoom cells on this level
      thisLevelGrid->unflagAllCells();
      thisLevelGrid->flagCells(sourceArray);
      auto proxyGrid = thisLevelGrid->makeProxyGridToMatch(*aboveLevelGrid);

      std::vector<size_t> zoomParticleArrayLevelAbove;
      proxyGrid->getFlaggedCells(zoomParticleArrayLevelAbove);

      // zoomParticleArrayLevelAbove now has the minimal mask.
      // Expand the minimal mask by nPadCells in each direction
      aboveLevelGrid->unflagAllCells();
      aboveLevelGrid->flagCells(zoomParticleArrayLevelAbove);
      aboveLevelGrid->expandFlaggedRegion(nPadCells);
      zoomParticleArrayLevelAbove.clear();
      aboveLevelGrid->getFlaggedCells(zoomParticleArrayLevelAbove);
      this->zoomParticleArray[level - 1] = zoomParticleArrayLevelAbove;
    }
    updateParticleMapper();
  }

  //! Define cell currently pointed at by coordinates
  void setCentre(T xin, T yin, T zin) {
    x0 = xin;
    y0 = yin;
    z0 = zin;
  }

  //! Returns the coordinate of the centre of the lowest resolution grid
  Coordinate<T> getBoxCentre() {
    T boxsize = multiLevelContext.getGridForLevel(0).thisGridSize;
    return Coordinate<T>(boxsize / 2, boxsize / 2, boxsize / 2);
  }

  //! Calculates and prints chi^2 for the underlying field
  virtual void getFieldChi2() {
    initialiseRandomComponentIfUninitialised();
    T val = this->outputFields[0]->getChi2();
    size_t dof = this->multiLevelContext.getNumDof();
    logging::entry() << "Calculated chi^2 = " << std::setprecision(10) << val << " (dof = " << dof << ")" << std::endl;
  }

  //! Calculate unmodified property for the dark matter field:
  virtual void calculate(std::string name) {
    initialiseRandomComponentIfUninitialised();

    GridDataType val = modificationManager.calculateCurrentValueByName(name, this->variance_filterscale);
    logging::entry() << name << ": calculated value = " << val << endl;
  }

  //! Define a modification to be applied to the field
  /*!
   * \param name String with the name of the modification to be added.
   * @param type  Modification can be relative to existing value or absolute
   * @param target Absolute target or factor by which the existing will be multiplied
   */
  virtual void modify(string name, string type, float target) {
    initialiseRandomComponentIfUninitialised();
    // needed for relative modifications where a current value will be evaluated

    modificationManager.addModificationToList(name, type, target, this->initial_number_steps, this->precision,
                                              this->variance_filterscale);
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
    if (this->outputFormat == io::OutputFormat::unknown) {
      throw (std::runtime_error("No output format specified!"));
    }

    initialiseRandomComponentIfUninitialised();

    if (modificationManager.hasModifications())
      applyModifications();

    write();
  }


  //! Reverses the sign of the field.
  virtual void reverse() {
    initialiseRandomComponentIfUninitialised();
    outputFields[0]->reverse();
  }

  //! Reverses the sign of the low-k modes.
  virtual void reverseSmallK(T kmax) {

    initialiseRandomComponentIfUninitialised();

    T k2max = kmax * kmax;

    for (size_t level = 0; level < multiLevelContext.getNumLevels(); ++level) {
      auto &field = outputFields[0]->getFieldForLevel(level);
      field.toFourier();

      field.forEachFourierCell([k2max](std::complex<T> val, T kx, T ky, T kz) {
        T k2 = kx * kx + ky * ky + kz * kz;
        if (k2 < k2max && k2 != 0) {
          val = -val;
        }
        return val;
      });
    }
  }

};

#endif
