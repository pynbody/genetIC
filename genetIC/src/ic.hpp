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

  fields::OutputField<GridDataType> outputField;
  modifications::ModificationManager<GridDataType> modificationManager;
  fields::RandomFieldGenerator<GridDataType> randomFieldGenerator;

  cosmology::CAMB<GridDataType> spectrum;

  //! DM supersampling to perform on zoom grid
  int supersample;

  //! Subsampling on base grid
  int subsample;

  T xOffOutput, yOffOutput, zOffOutput;


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


  std::vector<size_t> flaggedParticles;
  std::vector<std::vector<size_t>> zoomParticleArray;


	//! Coordinates of the cell of current interest
  T x0, y0, z0;

  shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pMapper;
  shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pInputMapper;
  shared_ptr<multilevelcontext::MultiLevelContextInformation<GridDataType>> pInputMultiLevelContext;

  shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>> pParticleGenerator;

  using FieldType = fields::Field<GridDataType>;

  tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter;


public:
  ICGenerator(tools::ClassDispatch<ICGenerator<GridDataType>, void> &interpreter) :
    outputField(multiLevelContext),
    modificationManager(multiLevelContext, cosmology, outputField),
    randomFieldGenerator(outputField),
    pMapper(new particle::mapper::ParticleMapper<GridDataType>()),
    interpreter(interpreter) {
    pInputMapper = nullptr;
    pInputMultiLevelContext = nullptr;
    cosmology.hubble = 0.701;   // old default
    cosmology.OmegaBaryons0 = -1.0;
    cosmology.ns = 0.96;      // old default
    cosmology.TCMB = 2.725;
    haveInitialisedRandomComponent = false;
    supersample = 1;
    subsample = 1;
    xOffOutput = 0;
    yOffOutput = 0;
    zOffOutput = 0;
    exactPowerSpectrum = false;
    allowStrayParticles=false;
    pParticleGenerator = std::make_shared<particle::NullMultiLevelParticleGenerator<GridDataType>>();
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

	// TODO What is this offset ?
  void offsetOutput(T x, T y, T z) {
    xOffOutput = x;
    yOffOutput = y;
    zOffOutput = z;
    updateParticleMapper();
  }

  void setSigma8(T in) {
    cosmology.sigma8 = in;
  }

	//TODO What is this doing ? What is in ?
  void setSupersample(int in) {
    supersample = in;
    updateParticleMapper();
  }

	//TODO What is this doing ? What is in ?
  void setSubsample(int in) {
    subsample = in;
    updateParticleMapper();
  }

	//! //TODO Is z0 the redhsift where Zeldo ends or starts ?
  void setZ0(T in) {
    cosmology.redshift = in;
    cosmology.scalefactor = 1. / (cosmology.redshift + 1.);
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

    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));

    addLevelToContext(spectrum, boxSize, n);
    updateParticleMapper();

  }

  void setns(T in) {
    cosmology.ns = in;
  }

  /*! Define the zoomed grid encompassing all flagged cells/particles
   * \param zoomfac //TODO
   * \param n Number of cells in the zoom grid
   */
  void initZoomGrid(size_t zoomfac, size_t n) {
    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));

    if (multiLevelContext.getNumLevels() < 1)
      throw std::runtime_error("Cannot initialise a zoom grid before initialising the base grid");

    grids::Grid<T> &gridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);
    int nAbove = (int) gridAbove.size;

    storeCurrentCellFlagsAsZoomMask(multiLevelContext.getNumLevels());
    vector<size_t> &newLevelZoomParticleArray = zoomParticleArray.back();

    if(newLevelZoomParticleArray.size()==0)
      throw std::runtime_error("Cannot initialise zoom without marking particles to be replaced");

    // find boundaries
    Window<int> zoomWindow(gridAbove.getEffectiveSimulationSize(), gridAbove.getCellCoordinate(newLevelZoomParticleArray[0]));

    for (auto cell_id : newLevelZoomParticleArray) {
      zoomWindow.expandToInclude(gridAbove.getCellCoordinate(cell_id));
    }

    int n_required = zoomWindow.getMaximumDimension();

    // Now see if the zoom the user chose is OK
    int n_user = nAbove / zoomfac;
    if (n_required>n_user && !allowStrayParticles) {
      throw (std::runtime_error(
        "Zoom particles do not fit in specified sub-box. Decrease zoom, or choose different particles"));
    }

    if(n_required<n_user && !allowStrayParticles)
      zoomWindow.expandSymmetricallyToSize(n_user);

    auto lci = zoomWindow.getLowerCornerInclusive();

    initZoomGridWithOriginAt(lci.x, lci.y, lci.z, zoomfac, n);

  }

  void storeCurrentCellFlagsAsZoomMask(size_t level) {
    assert(level>0);

    if(zoomParticleArray.size()<level)
      zoomParticleArray.emplace_back();

    assert(zoomParticleArray.size()>=level);

    grids::Grid<T> &gridAbove = multiLevelContext.getGridForLevel(level-1);

    vector<size_t> &levelZoomParticleArray = zoomParticleArray[level-1];
    levelZoomParticleArray.clear();
    gridAbove.getFlaggedCells(levelZoomParticleArray);
  }

	/*! Define a zoomed grid with user defined origin
	 * \param x0, y0, z0  //TODO Coordinates of centre of grid or upper left corner ?
	 * \param zoomfac //TODO What?
	 * \param n Number of cells in the zoom grid
	 */
  void initZoomGridWithOriginAt(int x0, int y0, int z0, size_t zoomfac, size_t n) {
    grids::Grid<T> &gridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels() - 1);
    int nAbove = int(gridAbove.size);

    storeCurrentCellFlagsAsZoomMask(multiLevelContext.getNumLevels());

    vector<size_t> trimmedParticleArray;
    vector<size_t> &newLevelZoomParticleArray = zoomParticleArray.back();
    int nCoarseCellsOfZoomGrid = nAbove / int(zoomfac);

    Coordinate<int> lowerCorner(x0,y0,z0);
    Coordinate<int> upperCornerExclusive = lowerCorner+nCoarseCellsOfZoomGrid;

    if(gridAbove.coversFullSimulation())
      upperCornerExclusive = gridAbove.wrapCoordinate(upperCornerExclusive);
    else {
      if (upperCornerExclusive.x>nAbove || upperCornerExclusive.y>nAbove || upperCornerExclusive.z>nAbove)
        throw std::runtime_error("Attempting to initialise a zoom grid that falls outside the parent grid");
    }


    size_t missed_particle=0;

    Window<int> zoomWindow = Window<int>(gridAbove.getEffectiveSimulationSize(),
                                         lowerCorner, upperCornerExclusive);


    // Make list of the particles, excluding those that fall outside the new high-res box. Alternatively,
    // if allowStrayParticles is true, keep even those outside the high-res box but report the number
    // in this category.
    for (size_t i = 0; i < newLevelZoomParticleArray.size(); i++) {
      bool include=true;
      auto coord = gridAbove.getCellCoordinate(newLevelZoomParticleArray[i]);
      if(!zoomWindow.contains(coord)) {
        missed_particle+=1;
        include = false;
      }

      if(include || allowStrayParticles) {
        trimmedParticleArray.push_back(newLevelZoomParticleArray[i]);
      }
    }

    if(missed_particle>0) {
      cerr << "WARNING: the requested zoom particles do not all fit in the requested zoom window" << endl;
      if(allowStrayParticles) {
        cerr << "         of " << newLevelZoomParticleArray.size() << " particles, " << missed_particle << " will be interpolated from LR grid (stray particle mode)" << endl;
      } else {
        cerr << "         of " << newLevelZoomParticleArray.size() << " particles, " << missed_particle << " have been omitted" << endl;
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
    cout << "  Origin in parent grid = " << lowerCorner << endl;
    cout << "  Low-left corner       = " << newGrid.offsetLower.x << ", " << newGrid.offsetLower.y << ", "
         << newGrid.offsetLower.z << endl;
    cout << "  Num particles         = " << newLevelZoomParticleArray.size() << endl;

    updateParticleMapper();

    cout << "  Total particles = " << pMapper->size() << endl;
  }

  virtual void
  addLevelToContext(const cosmology::CAMB<GridDataType> &spectrum, T size, size_t nside, const Coordinate<T> &offset = {0, 0, 0}) {
    // This forwards to multiLevelContext but is required because it is overriden in DummyICGenerator,
    // which needs to ensure that grids are synchronised between two different contexts
    multiLevelContext.addLevel(spectrum, size, nside, offset);
  }


  void setSeed(int seed) {
    randomFieldGenerator.seed(seed);
  }

  void setSeedFourier(int seed) {
    randomFieldGenerator.seed(seed);
    randomFieldGenerator.setDrawInFourierSpace(true);
    randomFieldGenerator.setReverseRandomDrawOrder(false);
  }

	//TODO What is this Reverse ?
  void setSeedFourierReverseOrder(int seed) {
    randomFieldGenerator.seed(seed);
    randomFieldGenerator.setDrawInFourierSpace(true);
    randomFieldGenerator.setReverseRandomDrawOrder(true);
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
	 * \param format //TODO Which int is which ?
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
  virtual void zeroLevel(int level) {
    cerr << "*** Warning: your script calls zeroLevel(" << level << "). This is intended for testing purposes only!"
         << endl;

    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    auto &fieldData = outputField.getFieldForLevel(level).getDataVector();
    std::fill(fieldData.begin(), fieldData.end(), 0);
  }


  virtual void applyPowerSpec() {
    if (this->exactPowerSpectrum) {
      outputField.enforceExactPowerSpectrum();
    } else {
      outputField.applyPowerSpectrum();
    }
  }

  template<typename TField>
  void dumpGridData(int level, const TField &data) {
    grids::Grid<T> &levelGrid = multiLevelContext.getGridForLevel(level);

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

	// TODO not quire sure what the point of this method is ? Looks like a debugging tool
  virtual void saveTipsyArray(string fname) {
    io::tipsy::saveFieldTipsyArray(fname, *pMapper, *pParticleGenerator, outputField);
  }

	//! Dumps field at a given level in a file named grid-level
  virtual void dumpGrid(int level = 0) {
    outputField.toReal();
    dumpGridData(level, outputField.getFieldForLevel(level));
  }

	// TODO Is this used at all ?
  virtual void dumpGridFourier(int level = 0) {
    fields::Field<complex<T>, T> fieldToWrite = tools::numerics::fourier::getComplexFourierField(
      outputField.getFieldForLevel(level));
    dumpGridData(level, fieldToWrite);
  }

	//! Dumps power spectrum generated from the field and the theory at a given level in a .ps file
  virtual void dumpPS(int level = 0) {
    auto &field = outputField.getFieldForLevel(level);
    field.toFourier();
    cosmology::dumpPowerSpectrum(field,
                                 multiLevelContext.getCovariance(level),
                                 (getOutputPath() + "_" + ((char) (level + '0')) + ".ps").c_str());
  }


  virtual void initialiseParticleGenerator() {
    // in principle this could now be easily extended to slot in higher order PT or other
    // methods of generating the particles from the fields

    using GridLevelGeneratorType = particle::ZeldovichParticleGenerator<GridDataType>;

    pParticleGenerator = std::make_shared<
      particle::MultiLevelParticleGenerator<GridDataType, GridLevelGeneratorType>>(outputField, cosmology);

  }

	// TODO Not quite sure what this if for ?
  void setInputMapper(std::string fname) {
    DummyICGenerator<GridDataType> pseudoICs(this);
    auto dispatch = interpreter.specify_instance(pseudoICs);
    ifstream inf;
    inf.open(fname);


    if (!inf.is_open())
      throw std::runtime_error("Cannot open IC paramfile for relative_to command");
    cerr << "******** Running commands in" << fname << " to work out relationship ***********" << endl;

    tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(fname));

    dispatch.run_loop(inf);
    cerr << *(pseudoICs.pMapper) << endl;
    cerr << "******** Finished with" << fname << " ***********" << endl;
    pInputMapper = pseudoICs.pMapper;
    pInputMultiLevelContext = std::make_shared<multilevelcontext::MultiLevelContextInformation<GridDataType>>
      (pseudoICs.multiLevelContext);
  }

  /*! Get the grid on which the output is defined for a particular level.
   *
   * This may differ from the grid on which the fields are defined either because there is an offset or
   * there are differences in the resolution between the output and the literal fields.
   *
   */
  std::shared_ptr<grids::Grid<T>> getOutputGrid(int level = 0) {
    auto gridForOutput = multiLevelContext.getGridForLevel(level).shared_from_this();

    if (xOffOutput != 0 || yOffOutput != 0 || zOffOutput != 0) {
      gridForOutput = std::make_shared<grids::OffsetGrid<T>>(gridForOutput,
                                                             xOffOutput, yOffOutput, zOffOutput);
    }
    if(allowStrayParticles && level>0) {
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
      // Grafic format just writes out the grids in turn
      pMapper = std::make_shared<particle::mapper::GraficMapper<GridDataType>>(multiLevelContext);
      return;
    }

    // make a basic mapper for the coarsest grid
    pMapper = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
      new particle::mapper::OneLevelParticleMapper<GridDataType>(
          getOutputGrid(0)
      ));


    if (nLevels >= 2) {

      for (size_t level = 1; level < nLevels; level++) {

        auto pFine = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
          new particle::mapper::OneLevelParticleMapper<GridDataType>(getOutputGrid(level)));

        pMapper = std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>>(
          new particle::mapper::TwoLevelParticleMapper<GridDataType>(pMapper, pFine, zoomParticleArray[level - 1]));
      }
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


  virtual void write() {
    using namespace io;

    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    initialiseParticleGenerator();

    cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
    cerr << (*pMapper);

    T boxlen = multiLevelContext.getGridForLevel(0).periodicDomainSize;

    switch (outputFormat) {
      case OutputFormat::gadget2:
      case OutputFormat::gadget3:
        gadget::save(getOutputPath() + ".gadget", boxlen, *pMapper,
                     *pParticleGenerator,
                     cosmology, static_cast<int>(outputFormat));
        break;
      case OutputFormat::tipsy:
        tipsy::save(getOutputPath() + ".tipsy", boxlen, *pParticleGenerator,
                    pMapper, cosmology);
        break;
      case OutputFormat::grafic:
        grafic::save(getOutputPath() + ".grafic", *pParticleGenerator, multiLevelContext, cosmology);
        break;
      default:
        throw std::runtime_error("Unknown output format");
    }

  }

  void initialiseRandomComponent() {
    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to re-draw the random field after it was already initialised"));

    randomFieldGenerator.draw();
    applyPowerSpec();

    haveInitialisedRandomComponent = true;
  }


protected:

  int deepestLevelWithParticlesSelected() {
    for (size_t i = multiLevelContext.getNumLevels() - 1; i >= 0; --i) {
      if (multiLevelContext.getGridForLevel(i).hasFlaggedCells())
        return i;
    }
    throw std::runtime_error("No level has any particles selected");
  }

  int deepestLevel() {
    //TODO: can this be removed?
    return multiLevelContext.getNumLevels();
  }

  T get_wrapped_delta(T x0, T x1) {
    return multiLevelContext.getGridForLevel(0).getWrappedOffset(x0, x1);
  }


  void getCentre() {
    x0 = 0;
    y0 = 0;
    z0 = 0;

    int level = deepestLevelWithParticlesSelected();

    std::vector<size_t> particleArray;
    grids::Grid<T> &grid = multiLevelContext.getGridForLevel(level);
    grid.getFlaggedCells(particleArray);

    auto p0_location = grid.getCellCentroid(particleArray[0]);

    for (size_t i = 0; i < particleArray.size(); i++) {
      auto pi_location = grid.getCellCentroid(particleArray[i]);
      x0 += get_wrapped_delta(pi_location.x, p0_location.x);
      y0 += get_wrapped_delta(pi_location.y, p0_location.y);
      z0 += get_wrapped_delta(pi_location.z, p0_location.z);
    }
    x0 /= particleArray.size();
    y0 /= particleArray.size();
    z0 /= particleArray.size();

    cerr << "Centre of region is " << setprecision(12) << x0 << " " << y0 << " " << z0 << endl;
  }


  void appendParticleIdFile(std::string filename) {

    cerr << "Loading " << filename << endl;

    io::getBuffer(flaggedParticles, filename);
    size_t size = flaggedParticles.size();
    std::sort(flaggedParticles.begin(), flaggedParticles.end());
    flaggedParticles.erase(std::unique(flaggedParticles.begin(), flaggedParticles.end()),
                           flaggedParticles.end());
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


	//TODO Is this loading existing partciles ?
  void loadID(string fname) {
    loadParticleIdFile(fname);
    getCentre();
  }

  void appendID(string fname) {
    appendParticleIdFile(fname);
    getCentre();
  }

  virtual void dumpID(string fname) {
    std::vector<size_t> results;
    cerr << "dumpID using current mapper:" << endl;
    cerr << (*pMapper);
    pMapper->getFlaggedParticles(results);
    io::dumpBuffer(results, fname);
  }

	//! TODO Is this finding the cell for a given particle ?
  void centreParticle(long id) {
    std::tie(x0, y0, z0) = multiLevelContext.getGridForLevel(0).getCellCentroid(id);
  }

	//! Flag the nearest cell to the coordinates currently pointed at
  void selectNearest() {
    auto &grid = multiLevelContext.getGridForLevel(deepestLevel() - 1);
    pMapper->unflagAllParticles();
    size_t id = grid.getCellContainingPoint(Coordinate<T>(x0, y0, z0));
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
        std::tie(xp, yp, zp) = grid->getCellCentroid(i);
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
	// TODO Is radius in pixel or Mpc ? Same for Cube
  void selectSphere(float radius) {
    T r2 = radius * radius;
    select([r2](T delta_x, T delta_y, T delta_z) -> bool {
      T r2_i = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
      return r2_i < r2;
    });

  }

	//! Flag all cells contained in the cube centered at the coordinates currently pointed at
  void selectCube(float side) {
    T side_by_2 = side / 2;
    select([side_by_2](T delta_x, T delta_y, T delta_z) -> bool {
      return abs(delta_x) < side_by_2 && abs(delta_y) < side_by_2 && abs(delta_z) < side_by_2;
    });
  }

	//! Define cell currently pointed at by coordinates
  void setCentre(T xin, T yin, T zin) {
    x0 = xin;
    y0 = yin;
    z0 = zin;
  }

  void calculate(string name) {
    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

		GridDataType val = modificationManager.calculateCurrentValueByName(name);

    cout << name << ": calculated value = " << val << endl;
  }

  void calculateVariance(T filterscale){
    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    GridDataType val = modificationManager.calculateVariance(filterscale);
    cout << "variance" << ": calculated value = " << val << endl;
  }

  virtual void modify(string name, string type, float target) {
    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

		modificationManager.addModificationToList(name, type, target);

  }

	virtual void quadraticallyModify(string name, string type, T target, int initNumberSteps, T precision, T filterscale){
		modificationManager.addQuadModificationToList(name, type, target, initNumberSteps, precision, filterscale);
	}

  void cov() {
  	modificationManager.print_covariance();
  }

	void clearModifications() {
		modificationManager.clearModifications();
	}


  virtual void applyModifications() {
    if (!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    modificationManager.applyModifications();
  }

  virtual void done() {
    T pre_modif_chi2 = outputField.getChi2();
    cerr << "BEFORE modifications chi^2=" << pre_modif_chi2 << endl;
    applyModifications();
    T post_modif_chi2 = outputField.getChi2();
    cerr << "AFTER  modifications chi^2=" << post_modif_chi2 << endl;
    cerr << "             delta-chi^2=" << post_modif_chi2 - pre_modif_chi2 << endl;
    write();
  }

	//TODO This method does not really belong here but in MultiLevelField class
  void reverse() {
    for_each_level(level) {
      auto &field = outputField.getFieldForLevel(level);
      size_t N = field.getGrid().size3;
      auto &field_data = field.getDataVector();
      for (size_t i = 0; i < N; i++)
        field_data[i] = -field_data[i];
    }
  }

	//TODO What is this for and why is it never tested ? Looks like inverted initial conditions properties
  void reseedSmallK(T kmax, int seed) {

    T k2max = kmax * kmax;


    // take a copy of all the fields
    std::vector<FieldType> fieldCopies;
    for_each_level(level) {
      auto &field = outputField.getFieldForLevel(level);
      field.toFourier();
      fieldCopies.emplace_back(field);
    }

    // remake the fields with the new seed
    randomFieldGenerator.seed(seed);
    initialiseRandomComponent();

    // copy back the old field
    for_each_level(level) {
      FieldType &oldField = fieldCopies[level];
      auto &field = outputField.getFieldForLevel(level);
      int k2max_i = tools::getRatioAndAssertInteger(k2max,field.getGrid().getFourierKmin());
      field.forEachFourierCell([k2max_i, &oldField](std::complex<T> val, int kx, int ky, int kz) {
        int k2_i = kx * kx + ky * ky + kz * kz;
        if (k2_i < k2max_i && k2_i != 0) {
          val = oldField.getFourierCoefficient(kx,ky,kz);
        }
        return val;
      });
    }

  }

  void reverseSmallK(T kmax) {

    T k2max = kmax * kmax;

    for_each_level(level) {
      auto &field = outputField.getFieldForLevel(level);
      field.toFourier();

      field.forEachFourierCell([k2max](std::complex<T> val, T kx, T ky, T kz) {
        T k2 = kx*kx+ky*ky+kz*kz;
        if (k2 < k2max && k2 != 0) {
          val = -val;
        }
        return val;
      });

      cerr << "reverseSmallK: k reversal at " << sqrt(k2max) << endl;
    }

  }


};

#endif
