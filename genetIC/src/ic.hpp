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

#include "src/cosmology/camb.hpp"
#include "src/simulation/window.hpp"

//TODO: remove ugly macro
#define for_each_level(level) for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level)

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

  //Vector of output fields (defaults to just a single field)
  std::vector<fields::OutputField<GridDataType>> outputField;
  int nFields;
  //fields::OutputField<GridDataType> outputField
  //modifications::ModificationManager<GridDataType> modificationManager;
  std::vector<modifications::ModificationManager<GridDataType>> modificationManager;
  //fields::RandomFieldGenerator<GridDataType> randomFieldGenerator;
  std::vector<fields::RandomFieldGenerator<GridDataType>> randomFieldGenerator;

  //Switch to re-direct calls to the baryon field to the DM field if we are not using it:
  std::vector<size_t> transferSwitch;



  cosmology::CAMB<GridDataType> spectrum;

  //! Gadget particle types to be generated on each level (default 1)
  std::vector<unsigned int> gadgetTypesForLevels;

  //! If true, at output time assign a different gadget particle type to the flagged particles on the finest known level
  bool flaggedParticlesHaveDifferentGadgetType;
  unsigned int flaggedGadgetParticleType;

  //! DM supersampling to perform on zoom grid
  int supersample;

  //! Subsampling on base grid
  int subsample;


  io::OutputFormat outputFormat;
  string outputFolder, outputFilename;

  //! Track whether the random realisation has yet been made
  //bool haveInitialisedRandomComponent;
  std::vector<bool> haveInitialisedRandomComponent;

  //! Enforce the exact power spectrum, as in Angulo & Pontzen 2016
  bool exactPowerSpectrum;

  /*!
   Stray" particles are high-res particles outside a high-res grid,
   constructed through interpolation of the surrounding low-res grid. Disabled by default.
   */
  bool allowStrayParticles;

  //! If true, the box are recentered on the last centered point in the parameter file
  //TODO This is currently a grafic only option and ought to be generic
  bool centerOnTargetRegion = false ;


  std::vector<size_t> flaggedParticles;
  std::vector<std::vector<size_t>> zoomParticleArray;


  //! Coordinates of the cell of current interest
  T x0, y0, z0;

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

  //shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>> pParticleGenerator;
  std::vector<shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>>> pParticleGenerator;

  using FieldType = fields::Field<GridDataType>;

  tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter;


public:
  ICGenerator(tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter) :
      //outputField(multiLevelContext),//This will have to be changed, because we can no longer provide it as an argument to modificationManager and randomFieldGenerator
      //modificationManager(multiLevelContext, cosmology, outputField),
      //randomFieldGenerator(outputField),
      pMapper(new particle::mapper::ParticleMapper<GridDataType>()),
      interpreter(interpreter) {

    //Set default parameters, in case input file didn't specify these:
    //By default, we assume there is only one field - the DM field:
    nFields = 1;
    outputField.emplace_back(multiLevelContext);
    modificationManager.emplace_back(multiLevelContext, cosmology, outputField[0]);
    randomFieldGenerator.emplace_back(outputField[0]);
    transferSwitch.push_back(0);

    pInputMapper = nullptr;
    pInputMultiLevelContext = nullptr;
    cosmology.hubble = 0.701;   // old default
    cosmology.OmegaBaryons0 = -1.0;
    cosmology.ns = 0.96;      // old default
    cosmology.TCMB = 2.725;
    //haveInitialisedRandomComponent = false;
    haveInitialisedRandomComponent.push_back(false);
    supersample = 1;
    subsample = 1;
    exactPowerSpectrum = false;
    allowStrayParticles = false;
    auto _pParticleGenerator = std::make_shared<particle::NullMultiLevelParticleGenerator<GridDataType>>();
    pParticleGenerator.push_back(_pParticleGenerator);
    flaggedParticlesHaveDifferentGadgetType = false;
    flaggedGadgetParticleType = 1;

  }

  ~ICGenerator() {

  }

  void setOmegaM0(T in) {
    cosmology.OmegaM0 = in;
  }

  void setTCMB(T in) {
    cosmology.TCMB = in;
  }

  void setOmegaB0(T in) {
    cosmology.OmegaBaryons0 = in;
    // now that we have gas, mapper may have changed:
    updateParticleMapper();
  }

  void setOmegaLambda0(T in) {
    cosmology.OmegaLambda0 = in;
  }

  void setHubble(T in) {
    cosmology.hubble = in;
  }

  //! Enables the use of stray particles
  void setStraysOn() {
    allowStrayParticles = true;
  }

  //! Enables the use of a baryon density field:
  void setUsingBaryons()
  {
    //Really have to make sure that this is run near the beginning of the parameter file,
    //otherwise we will encounter serious problems if we start doing calculations and then try to add in baryons.
    if(this->nFields > 1)
    {
        std::cerr<< "Already using baryons!" << std::endl;
    }
    else
    {
        this->nFields = 2;
        std::cerr<< "Baryon field enabled." << std::endl;
        outputField.emplace_back(multiLevelContext);
        modificationManager.emplace_back(multiLevelContext, cosmology, outputField[1]);
        randomFieldGenerator.emplace_back(outputField[1]);
        transferSwitch.push_back(1);
        haveInitialisedRandomComponent.push_back(false);
        auto _pParticleGenerator = std::make_shared<particle::NullMultiLevelParticleGenerator<GridDataType>>();
        pParticleGenerator.push_back(_pParticleGenerator);
        this->spectrum.enableAllTransfers();
    }
  }


  void setSigma8(T in) {
    cosmology.sigma8 = in;
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

  void setCenteringOnRegion(){
    this->centerOnTargetRegion = true;
    updateParticleMapper();
  }

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

    for(size_t i = 0;i < haveInitialisedRandomComponent.size();i++)
    {
        if(haveInitialisedRandomComponent[i])
        {
            throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));
        }
    }


    addLevelToContext(spectrum, boxSize, n);
    updateParticleMapper();

  }

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
    for(size_t i = 0;i < haveInitialisedRandomComponent.size();i++)
    {
        if(haveInitialisedRandomComponent[i])
        {
            throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));
        }
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

    if (n_required < n_user && !allowStrayParticles)
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
  }

  virtual void
  addLevelToContext(const cosmology::CAMB<GridDataType> &spectrum, T size, size_t nside,
                    const Coordinate<T> &offset = {0, 0, 0}) {
    // This forwards to multiLevelContext but is required because it is overriden in DummyICGenerator,
    // which needs to ensure that grids are synchronised between two different contexts
    multiLevelContext.addLevel(spectrum, size, nside, offset);
    this->gadgetTypesForLevels.push_back(1);
  }

  //!Set all fields to be randomised with the same seed.
  void setSeedAll(int seed)
  {
    for(size_t i = 0;i < outputField.size();i++)
    {
        this->setSeed(seed,i);
    }
  }
  //!Seed all fields in Fourier space with the same seed.
  void setSeedFourierAll(int seed)
  {
    for(size_t i = 0;i < outputField.size();i++)
    {
        this->setSeedFourier(seed,i);
    }
  }

  //Could do this with a default parameter too, but that plays
  //havoc with the recursive templates used by the parser,
  //which can't see the default arguments and thinks
  //not enough arguments have been supplied.



  //!Set the seed in real space.
  void setSeed(int seed,size_t nField) {
    randomFieldGenerator[nField].seed(seed);
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  void setSeed(int seed)
  {
    this->setSeed(seed,0);
  }



  //!Set the seed in fourier space
  void setSeedFourier(int seed,size_t nField) {
    randomFieldGenerator[nField].seed(seed);
    randomFieldGenerator[nField].setDrawInFourierSpace(true);
    randomFieldGenerator[nField].setReverseRandomDrawOrder(false);
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  void setSeedFourier(int seed)
  {
    this->setSeedFourier(seed,0);
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

  //! Reverses the order of draws between real and imaginary part of complex numbers
  /*!
   * Provided for compatibility problems as different compilers handle the draw order differently
   */
  void setSeedFourierReverseOrder(int seed,size_t nField) {
    randomFieldGenerator[nField].seed(seed);
    randomFieldGenerator[nField].setDrawInFourierSpace(true);
    randomFieldGenerator[nField].setReverseRandomDrawOrder(true);
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  void setSeedFourierReverseOrder(int seed)
  {
    this->setSeedFourierReverseOrder(seed,0);
  }

  //!Seed in reverse order for all fields with the same seed.
  void setSeedFourierReverseOrderAll(int seed)
  {
    for(size_t i = 0;i < outputField.size();i++)
    {
        this->setSeedFourierReverseOrder(seed,i);
    }
  }

  void setExactPowerSpectrumEnforcement() {
    exactPowerSpectrum = true;
  }

  //! Obtain power spectrum from a CAMB data file
  void setCambDat(std::string cambFilePath) {
    spectrum.read(cambFilePath, cosmology);
  }

  void setOutDir(std::string outputPath) {
    outputFolder = outputPath;
  }

  void setOutName(std::string outputFilename_) {
    outputFilename = outputFilename_;
  }

  //! Set formats from handled formats in io namespace
  /*!
   * \param format  2 =  Gadget2, 3 = Gadget3, 4 = tipsy, 5 = grafic
   */
  void setOutputFormat(int format) {
    outputFormat = static_cast<io::OutputFormat>(format);
    updateParticleMapper();
  }

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

  //! Zeroes field values at a given level. Meant for debugging.
  virtual void zeroLevel(size_t level,size_t nField) {
    cerr << "*** Warning: your script calls zeroLevel(" << level << "). This is intended for testing purposes only!"
         << endl;

    if (!haveInitialisedRandomComponent[nField])
      initialiseRandomComponent(nField);

    auto &fieldData = outputField[nField].getFieldForLevel(level).getDataVector();
    std::fill(fieldData.begin(), fieldData.end(), 0);
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  virtual void zeroLevel(size_t level)
  {
    this->zeroLevel(level,0);
  }

  //!Imports random field data for a level from a supplies file.
  virtual void importLevel(size_t level, std::string filename,size_t nField = 0) {
    if (!haveInitialisedRandomComponent[nField])
      initialiseRandomComponent(nField);

    cerr << "Importing random field on level " << level << " from " <<filename << endl;

    auto & levelField = outputField[nField].getFieldForLevel(level);
    levelField.loadGridData(filename);
    levelField.setFourier(false);
    cerr << "... success!" << endl;
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  virtual void importLevel(size_t level,std::string filename)
  {
    this->importLevel(level,filename,0);
  }

  //!Apply the power spectrum to the white noise fields
  virtual void applyPowerSpec(size_t nField) {
    if (this->exactPowerSpectrum) {
      outputField[nField].enforceExactPowerSpectrum();
    } else {
      outputField[nField].applyPowerSpectrum();
    }
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  virtual void applyPowerSpec()
  {
    this->applyPowerSpec(0);
  }

  template<typename TField>
  void dumpGridData(size_t level, const TField &data) {


  //Unsure if this is necessary. May have to do this only for a single field,
  //But as it is supplied as an argument, we can't tell which one it was!
    for(size_t i = 0;i < modificationManager.size();i++)
    {
        if (!haveInitialisedRandomComponent[i])
            initialiseRandomComponent(i);
    }

    auto levelGrid = data.getGrid();

    ostringstream filename;
    filename << outputFolder << "/grid-" << level << ".npy";

    data.dumpGridData(filename.str());

    filename.str("");

    filename << outputFolder << "/grid-info-" << level << ".txt";

    ofstream ifile;
    ifile.open(filename.str());
    cerr << "Writing to " << filename.str() << endl;

    ifile << levelGrid.offsetLower.x << " " << levelGrid.offsetLower.y << " "
          << levelGrid.offsetLower.z << " " << levelGrid.thisGridSize << endl;
    ifile << "The line above contains information about grid level " << level << endl;
    ifile << "It gives the x-offset, y-offset and z-offset of the low-left corner and also the box length" << endl;
    ifile.close();
  }

  //! Dumps field in a tipsy format
  virtual void saveTipsyArray(string fname,size_t nField = 0) {
    io::tipsy::saveFieldTipsyArray(fname, *pMapper, *pParticleGenerator[nField], outputField[nField]);
  }

  //! Dumps field at a given level in a file named grid-level
  virtual void dumpGrid(size_t level,size_t nField) {
    //May need to modify this to work better with multiple fields - do we want them in the same or different fields?
    if(level < 0 || level > this->multiLevelContext.getNumLevels() -1){
      throw std::runtime_error("Trying to dump an undefined level");
    }
    outputField[nField].toReal();
    dumpGridData(level, outputField[nField].getFieldForLevel(level));
    outputField[nField].toFourier();
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  virtual void dumpGrid(size_t level) {
    this->dumpGrid(level,0);
  }

  // TODO Is this used at all ? Should be linked to a command in main if we want to keep it.
  virtual void dumpGridFourier(size_t level = 0,size_t nField = 0) {
    fields::Field<complex<T>, T> fieldToWrite = tools::numerics::fourier::getComplexFourierField(
        outputField[nField].getFieldForLevel(level));
    dumpGridData(level, fieldToWrite);
  }

  //! Dumps power spectrum generated from the field and the theory at a given level in a .ps file
  virtual void dumpPS(size_t level,size_t nField) {
    auto &field = outputField[nField].getFieldForLevel(level);
    field.toFourier();
    cosmology::dumpPowerSpectrum(field,
                                 multiLevelContext.getCovariance(level,nField),
                                 (getOutputPath() + "_" + ((char) (level + '0')) + ".ps").c_str());
  }

  //!For backwards compatibility with scripts that don't use baryon transfer function
  virtual void dumpPS(size_t level) {
    this->dumpPS(level,0);
  }
  //!For backwards compatibility with scripts that don't use baryon transfer function
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
      dumpGridData(level, maskfield->getFieldForLevel(level));
    }
  }


  virtual void initialiseParticleGenerator(size_t nField = 0) {
    // in principle this could now be easily extended to slot in higher order PT or other
    // methods of generating the particles from the fields

    using GridLevelGeneratorType = particle::ZeldovichParticleGenerator<GridDataType>;

    pParticleGenerator[nField] = std::make_shared<
        particle::MultiLevelParticleGenerator<GridDataType, GridLevelGeneratorType>>(outputField[nField], cosmology);

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
      auto gasMapper = pMapper->addGas(cosmology.OmegaBaryons0 / cosmology.OmegaM0,
                                       {multiLevelContext.getGridForLevel(nLevels - 1).shared_from_this()});

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

  void reflag() {

    if (pInputMapper != nullptr) {
      pMapper->unflagAllParticles();
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
    //Some work needed here to make sure this works correctly with multiple fields.

    for(size_t i = 0;i < haveInitialisedRandomComponent.size();i++)
    {
        if (!haveInitialisedRandomComponent[i])
            initialiseRandomComponent(i);

        initialiseParticleGenerator(i);
    }


    cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
    cerr << (*pMapper);

    T boxlen = multiLevelContext.getGridForLevel(0).periodicDomainSize;

    switch (outputFormat) {
      case OutputFormat::gadget2:
      case OutputFormat::gadget3:
        gadget::save(getOutputPath() + ".gadget", boxlen, *pMapper,
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
                       subsample, supersample, zoomParticleArray,outputField);
        } else {
          grafic::save(getOutputPath() + ".grafic",
                       pParticleGenerator, multiLevelContext, cosmology, pvarValue, this->getBoxCentre(),
                       subsample, supersample, zoomParticleArray,outputField);
        }

        break;
      default:
        throw std::runtime_error("Unknown output format");
    }

  }

  void initialiseRandomComponent(size_t nField = 0) {
    if (haveInitialisedRandomComponent[nField])
      throw (std::runtime_error("Trying to re-draw the random field after it was already initialised"));

    randomFieldGenerator[nField].draw();
    applyPowerSpec(nField);

    haveInitialisedRandomComponent[nField] = true;
  }


protected:

  size_t deepestLevelWithParticlesSelected() {
    return multiLevelContext.deepestLevelwithFlaggedCells();
  }

  size_t deepestLevel() {
    return multiLevelContext.getNumLevels();
  }

  T get_wrapped_delta(T x0, T x1) {
    return multiLevelContext.getGridForLevel(0).getWrappedOffset(x0, x1);
  }


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
  }

  //! Flag the nearest cell to the coordinates currently pointed at
  void selectNearest() {
    auto &grid = multiLevelContext.getGridForLevel(deepestLevel() - 1);
    pMapper->unflagAllParticles();
    size_t id = grid.getIndexFromPoint(Coordinate<T>(x0, y0, z0));
    cerr << "selectNearest " << x0 << " " << y0 << " " << z0 << " " << id << " " << endl;
    grid.flagCells({id});

  }

  void select(std::function<bool(T, T, T)> inclusionFunction) {
    T delta_x, delta_y, delta_z;
    T xp, yp, zp;

    flaggedParticles.clear();


    // unflag all grids first. This can't be in the loop below in case there are subtle
    // relationships between grids (in particular the ResolutionMatchingGrid which actually
    // points to two levels simultaneously).
    for_each_level(level) {
      getOutputGrid(level)->unflagAllCells();
    }

    for_each_level(level) {
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
    for_each_level(level) {
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
  }

  Coordinate<T> getBoxCentre(){
    T boxsize = multiLevelContext.getGridForLevel(0).thisGridSize;
    return Coordinate<T>(boxsize/2,boxsize/2,boxsize/2);
  }

  virtual void getFieldChi2(size_t nField = 0){
    if (!haveInitialisedRandomComponent[nField])
      initialiseRandomComponent(nField);

    std::cerr << "Calculated chi^2 = " << this->outputField[nField].getChi2() <<std::endl;
  }

  //! Calculate physical quantities of the field
  /*!
   * @param filterscale Filtering scale in Mpc if the quantity needs it
   */
  void calculate(string name,size_t nField = 0) {
    if (!haveInitialisedRandomComponent[nField])
      initialiseRandomComponent(nField);

    GridDataType val = modificationManager[nField].calculateCurrentValueByName(name, this->variance_filterscale);

    cout << name << ": calculated value = " << val << endl;
  }

  //! Define a modification to be applied to the field
  /*!
   * @param type  Modification can be relative to existing value or absolute
   * @param target Absolute target or factor by which the existing will be multiplied
   */
  virtual void modify(string name, string type, float target,size_t nField = 0) {

    if (!haveInitialisedRandomComponent[nField])
      initialiseRandomComponent(nField);

    modificationManager[nField].addModificationToList(name, type, target,this->initial_number_steps, this->precision, this->variance_filterscale);
  }

  //! Empty modification list
  void clearModifications() {
    for(size_t i = 0;i < modificationManager.size();i++)
    {
        this->clearModifications(i);
    }
  }

  void clearModifications(size_t nField)
  {
    modificationManager[nField].clearModifications();
  }


  //! Apply the algorithm to produce the modified field
  virtual void applyModifications() {
    //modificationManager.applyModifications();
    for(size_t i = 0;i < modificationManager.size();i++)
    {
        this->applyModifications(i);
    }
  }

  virtual void applyModifications(size_t nField)
  {
    modificationManager[nField].applyModifications();
  }

  //! Apply the modifications, calculate the corresponding delta chi^2 and recombine low and high-ks between grids to write a particle output
  virtual void done() {
  /*
    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    if(modificationManager.hasModifications())
      applyModifications();

    write();
    */
    for(size_t i = 0;i < modificationManager.size();i++)
    {
        this->done(i);
    }

    //Not sure how well this works if we call it multiple times yet - need to test.
    //Just call it on the DM particles for now.
    write();
  }

  virtual void done(size_t nField)
  {
    if (!haveInitialisedRandomComponent[nField])
      initialiseRandomComponent(nField);

    if(modificationManager[nField].hasModifications())
      applyModifications(nField);

      //May need to edit this in the case of multiple fields!
    //write(nField);
  }

  //TODO This method does not really belong here but in MultiLevelField class
  void reverse(size_t nField = 0) {
    for_each_level(level) {
      auto &field = outputField[nField].getFieldForLevel(level);
      size_t N = field.getGrid().size3;
      auto &field_data = field.getDataVector();
      for (size_t i = 0; i < N; i++)
        field_data[i] = -field_data[i];
    }
  }

  //TODO What is this method doing? Looks like inverted initial conditions properties
  void reseedSmallK(T kmax, int seed,size_t nField = 0) {

    T k2max = kmax * kmax;


    // take a copy of all the fields
    std::vector<FieldType> fieldCopies;
    for_each_level(level) {
      auto &field = outputField[nField].getFieldForLevel(level);
      field.toFourier();
      fieldCopies.emplace_back(field);
    }

    // remake the fields with the new seed
    randomFieldGenerator[nField].seed(seed);
    initialiseRandomComponent(nField);

    // copy back the old field
    for_each_level(level) {
      FieldType &oldField = fieldCopies[level];
      auto &field = outputField[nField].getFieldForLevel(level);
      int k2max_i = tools::getRatioAndAssertInteger(k2max, field.getGrid().getFourierKmin());
      field.forEachFourierCell([k2max_i, &oldField](std::complex<T> val, int kx, int ky, int kz) {
        int k2_i = kx * kx + ky * ky + kz * kz;
        if (k2_i < k2max_i && k2_i != 0) {
          val = oldField.getFourierCoefficient(kx, ky, kz);
        }
        return val;
      });
    }

  }

  void reverseSmallK(T kmax,size_t nField = 0) {

    T k2max = kmax * kmax;

    for_each_level(level) {
      auto &field = outputField[nField].getFieldForLevel(level);
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

  //Functions to control whether we use only the dark matter or include baryons:
  //Switch to using DM only:
  /*
  void setDMOnly()
  {
    this->spectrum.setDMOnly();
  }
  //Switch to allowing baryon transfer functions:
  void enableAllTransfers()
  {
    this->spectrum.enableAllTransfers();
  }
  */
};

#endif
