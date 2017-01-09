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

#include "numpy.hpp"
#include "onelevelconstraint.hpp"
#include "constraintapplicator.hpp"
#include "filesystem.h"
#include "randomfieldgenerator.hpp"
#include "MultiLevelConstraintGenerator.h"
#include "multilevelfield.hpp"
#include "particles/generator.hpp"
#include "particles/multilevelgenerator.hpp"
#include "parser.hpp"

//TODO: remove ugly macro
#define for_each_level(level) for(size_t level=0; level<multiLevelContext.getNumLevels(); ++level)

using namespace std;

template<typename T>
class DummyICGenerator;


template<typename T>
class ICGenerator {
  /** The top level object responsible for coordinating the generation of initial conditions, including GM constraints.
   *
   * This class exposes all the methods which are made scriptable by main.cpp
   */
protected:

  using GridPtrType = std::shared_ptr<Grid<T>>;

  friend class DummyICGenerator<T>;

  CosmologicalParameters<T> cosmology;
  MultiLevelContextInformation<T> multiLevelContext;
  OutputField<std::complex<T>> outputField;
  ConstraintApplicator<T> constraintApplicator;
  MultiLevelConstraintGenerator<T> constraintGenerator;
  RandomFieldGenerator<std::complex<T>> randomFieldGenerator;

  CAMB<T> spectrum;

  int supersample, subsample;               // DM supersampling to perform on zoom grid, and subsampling on base grid

  T xOffOutput, yOffOutput, zOffOutput;


  io::OutputFormat outputFormat;
  string outputFolder, outputFilename;

  bool haveInitialisedRandomComponent;
  bool exactPowerSpectrum;

  std::vector<size_t> flaggedParticles;
  std::vector<std::vector<size_t>> zoomParticleArray;


  T x0, y0, z0;

  shared_ptr<particle::ParticleMapper<T>> pMapper;
  shared_ptr<particle::ParticleMapper<T>> pInputMapper;

  shared_ptr<particle::AbstractMultiLevelParticleGenerator<T>> pParticleGenerator;

  using RefFieldType = std::vector<std::complex<T>> &;
  using FieldType = std::vector<std::complex<T>>;

  ClassDispatch<ICGenerator<T>, void> &interpreter;


public:
  ICGenerator(ClassDispatch<ICGenerator<T>, void> &interpreter) :
                                                      outputField(multiLevelContext),
                                                      constraintApplicator(&multiLevelContext, &outputField),
                                                      constraintGenerator(multiLevelContext, cosmology),
                                                      randomFieldGenerator(outputField),
                                                      pMapper(new particle::ParticleMapper<T>()),
                                                      interpreter(interpreter)
  {
    pInputMapper = nullptr;
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
    pParticleGenerator = std::make_shared<particle::NullMultiLevelParticleGenerator<T>>();
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

  void offsetOutput(T x, T y, T z) {
    xOffOutput = x;
    yOffOutput = y;
    zOffOutput = z;
    updateParticleMapper();
  }

  void setSigma8(T in) {
    cosmology.sigma8 = in;
  }

  void setSupersample(int in) {
    supersample = in;
    updateParticleMapper();
  }

  void setSubsample(int in) {
    subsample = in;
    updateParticleMapper();
  }

  void setZ0(T in) {
    cosmology.redshift = in;
    cosmology.scalefactor = 1. / (cosmology.redshift + 1.);
  }


  virtual void initBaseGrid(T boxSize, size_t n)
  {
    assert(boxSize>0);

    if(multiLevelContext.getNumLevels()>0)
      throw std::runtime_error("Cannot re-initialize the base grid");

    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));

    addLevelToContext(spectrum, boxSize, n);
    updateParticleMapper();

  }

  void setns(T in) {
    cosmology.ns = in;
  }

  void initZoomGrid(size_t zoomfac, size_t n) {
    if (haveInitialisedRandomComponent)
      throw (std::runtime_error("Trying to initialize a grid after the random field was already drawn"));

    if(multiLevelContext.getNumLevels()<1)
      throw std::runtime_error("Cannot initialise a zoom grid before initialising the base grid");

    Grid<T> & gridAbove = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels()-1);
    int nAbove = gridAbove.size;

    zoomParticleArray.emplace_back();
    vector<size_t> & newLevelZoomParticleArray = zoomParticleArray.back();
    gridAbove.getFlaggedCells(newLevelZoomParticleArray);

    // find boundaries
    int x0, x1, y0, y1, z0, z1;
    int x, y, z;

    x0 = y0 = z0 = gridAbove.size;
    x1 = y1 = z1 = 0;

    // TO DO: wrap the box sensibly

    for (size_t i = 0; i < newLevelZoomParticleArray.size(); i++) {
      std::tie(x, y, z) = gridAbove.getCellCoordinate(newLevelZoomParticleArray[i]);
      if (x < x0) x0 = x;
      if (y < y0) y0 = y;
      if (z < z0) z0 = z;
      if (x > x1) x1 = x;
      if (y > y1) y1 = y;
      if (z > z1) z1 = z;
    }

    // Now see if the zoom the user chose is OK
    int n_user = nAbove / zoomfac;
    if ((x1 - x0) > n_user || (y1 - y0) > n_user || (z1 - z0) > n_user) {
      throw (std::runtime_error(
        "Zoom particles do not fit in specified sub-box. Decrease zoom, or choose different particles. (NB wrapping not yet implemented)"));
    }
    // At this point we know things fit. All we need to do is choose
    // the correct offset to get the particles near the centre of the
    // zoom box.

    // Here is the bottom left of the box:
    x = (x0 + x1) / 2 - nAbove / (2 * zoomfac);
    y = (y0 + y1) / 2 - nAbove / (2 * zoomfac);
    z = (z0 + z1) / 2 - nAbove / (2 * zoomfac);

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (z < 0) z = 0;
    if (x > nAbove - nAbove / zoomfac) x = nAbove - nAbove / zoomfac;
    if (y > nAbove - nAbove / zoomfac) y = nAbove - nAbove / zoomfac;
    if (z > nAbove - nAbove / zoomfac) z = nAbove - nAbove / zoomfac;

    Coordinate<T> newOffsetLower = gridAbove.offsetLower + Coordinate<T>(x,y,z)*gridAbove.dx;

    addLevelToContext(spectrum, gridAbove.boxsize/zoomfac, n, newOffsetLower);

    Grid<T> &newGrid = multiLevelContext.getGridForLevel(multiLevelContext.getNumLevels()-1);

    cout << "Initialized a zoom region:" << endl;
    cout << "  Subbox length   = " << newGrid.boxsize << " Mpc/h" << endl;
    cout << "  n               = " << newGrid.size << endl;
    cout << "  dx              = " << newGrid.dx << endl;
    cout << "  Zoom factor     = " << zoomfac << endl;
    cout << "  Low-left corner = " << newGrid.offsetLower.x << ", " << newGrid.offsetLower.y << ", " << newGrid.offsetLower.z << endl;
    cout << "  Num particles   = " << newLevelZoomParticleArray.size() << endl;

    updateParticleMapper();

    cout << "  Total particles = " << pMapper->size() << endl;

  }

  virtual void addLevelToContext(const CAMB<T> &spectrum, T size, size_t nside, const Coordinate<T> & offset={0,0,0}) {
    // This forwards to multiLevelContext but is required because it is overriden in DummyICGenerator,
    // which needs to ensure that grids are synchronised between two different contexts
    multiLevelContext.addLevel(spectrum,size,nside,offset);
  }


  void setSeed(int in) {
    randomFieldGenerator.seed(in);
  }

  void setSeedFourier(int in) {
    randomFieldGenerator.seed(in);
    randomFieldGenerator.setDrawInFourierSpace(true);
    randomFieldGenerator.setReverseRandomDrawOrder(false);
  }

  void setSeedFourierReverseOrder(int in) {
    randomFieldGenerator.seed(in);
    randomFieldGenerator.setDrawInFourierSpace(true);
    randomFieldGenerator.setReverseRandomDrawOrder(true);
  }

  void setExactPowerSpectrumEnforcement() {
    exactPowerSpectrum = true;
  }

  void setCambDat(std::string in) {
    spectrum.read(in, cosmology);
  }

  void setOutDir(std::string in) {
    outputFolder = in;
  }

  void setOutName(std::string in) {
    outputFilename = in;
  }

  void setOutputFormat(int in) {
    outputFormat = static_cast<io::OutputFormat>(in);
    updateParticleMapper();
  }

  string getOutputPath() {
    ostringstream fname_stream;
    if (outputFilename.size() == 0) {
      fname_stream << outputFolder << "/IC_" << floatinfo<T>::name << "_z" << cosmology.redshift << "_" << multiLevelContext.getGridForLevel(0).size;
    } else {
      fname_stream << outputFolder << "/" << outputFilename;
    }
    return fname_stream.str();
  }

  virtual void zeroLevel(int level) {
    cerr << "*** Warning: your script calls zeroLevel(" << level << "). This is intended for testing purposes only!" << endl;

    if(!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    auto &fieldData = outputField.getFieldForLevel(level).getDataVector();
    std::fill(fieldData.begin(), fieldData.end(), 0);
  }


  virtual void applyPowerSpec() {
    if(this->exactPowerSpectrum) {
      outputField.enforceExactPowerSpectrum();
    } else {
      outputField.applyPowerSpectrum();
    }
  }

  template<typename TField>
  void dumpGridData(int level, const TField &data) {
    Grid<T> & levelGrid = multiLevelContext.getGridForLevel(level);
    assert(data.size() == levelGrid.size3);
    ostringstream filename;
    filename << outputFolder << "/grid-" << level << ".npy";

    int n = levelGrid.size;

    const int dim[3] = {n,n,n};
    numpy::SaveArrayAsNumpy(filename.str(), false, 3, dim, data.data());

    filename.str("");
    filename << outputFolder << "/grid-info-" << level << ".txt";

    ofstream ifile;
    ifile.open(filename.str());
    cerr << "Writing to " << filename.str() << endl;

    ifile << levelGrid.offsetLower.x << " " << levelGrid.offsetLower.y << " "
          << levelGrid.offsetLower.z << " " << levelGrid.boxsize << endl;
    ifile << "The line above contains information about grid level " << level << endl;
    ifile << "It gives the x-offset, y-offset and z-offset of the low-left corner and also the box length" << endl;
    ifile.close();
  }

  virtual void saveTipsyArray(string fname) {
    io::tipsy::saveFieldTipsyArray(fname, *pMapper, *pParticleGenerator, outputField);
  }

  virtual void dumpGrid(int level = 0) {
    outputField.toReal();
    dumpGridData(level, outputField.getFieldForLevel(level).getDataVector());
  }

  virtual void dumpPS(int level = 0) {
    auto & field = outputField.getFieldForLevel(level);
    field.toFourier();
    powsp_noJing(field,
                 multiLevelContext.getCovariance(level),
                 (getOutputPath() + "_" + ((char) (level + '0')) + ".ps").c_str(), field.getGrid().boxsize);
  }


  virtual void initialiseParticleGenerator() {
    // in principle this could now be easily extended to slot in higher order PT or other
    // methods of generating the particles from the fields

    using GridLevelGeneratorType = particle::ZeldovichParticleGenerator<T>;

    pParticleGenerator = std::make_shared<
      particle::MultiLevelParticleGenerator<T, GridLevelGeneratorType>>(outputField, cosmology);

  }

  void setInputMapper(std::string fname) {
    DummyICGenerator<T> pseudoICs(this);
    auto dispatch = interpreter.specify_instance(pseudoICs);
    ifstream inf;
    inf.open(fname);


    if (!inf.is_open())
      throw std::runtime_error("Cannot open IC paramfile for relative_to command");
    cerr << "******** Running commands in" << fname << " to work out relationship ***********" << endl;

    ChangeCwdWhileInScope temporary(getDirectoryName(fname));

    dispatch.run_loop(inf);
    cerr << *(pseudoICs.pMapper) << endl;
    cerr << "******** Finished with" << fname << " ***********" << endl;
    pInputMapper = pseudoICs.pMapper;

  }

  std::shared_ptr<Grid<T>> getGridWithOutputOffset(int level = 0) {
    auto gridForOutput = multiLevelContext.getGridForLevel(level).shared_from_this();

    if (xOffOutput != 0 || yOffOutput != 0 || zOffOutput != 0) {
      gridForOutput = std::make_shared<OffsetGrid<T>>(gridForOutput,
                                                      xOffOutput, yOffOutput, zOffOutput);
    }
    return gridForOutput;
  }

  void updateParticleMapper() {
    size_t nLevels = multiLevelContext.getNumLevels();

    if (nLevels == 0)
      return;

    // make a basic mapper for the coarsest grid
    pMapper = std::shared_ptr<particle::ParticleMapper<T>>(new particle::OneLevelParticleMapper<T>(
      getGridWithOutputOffset(0)
    ));


    if (nLevels >= 2) {

      for(size_t level=1; level<nLevels; level++) {

        auto pFine = std::shared_ptr<particle::ParticleMapper<T>>(
          new particle::OneLevelParticleMapper<T>(getGridWithOutputOffset(level)));

        pMapper = std::shared_ptr<particle::ParticleMapper<T>>(
          new particle::TwoLevelParticleMapper<T>(pMapper, pFine, zoomParticleArray[level-1]));
      }
    }

    if (cosmology.OmegaBaryons0 > 0) {

      // Add gas only to the deepest level. Pass the whole pGrid
      // vector if you want to add gas to every level.
      auto gasMapper = pMapper->addGas(cosmology.OmegaBaryons0 / cosmology.OmegaM0,
                                       {multiLevelContext.getGridForLevel(nLevels-1).shared_from_this()});

      bool gasFirst = outputFormat == io::OutputFormat::tipsy;

      // graft the gas particles onto the start of the map
      if (gasFirst)
        pMapper = std::make_shared<particle::AddGasMapper<T>>(
          gasMapper.first, gasMapper.second, true);
      else
        pMapper = std::make_shared<particle::AddGasMapper<T>>(
          gasMapper.second, gasMapper.first, false);

    }

    // potentially resample the lowest-level DM grid. Again, this is theoretically
    // more flexible if you pass in other grid pointers.
    if (supersample > 1)
      pMapper = pMapper->superOrSubSampleDM(supersample, {multiLevelContext.getGridForLevel(nLevels-1).shared_from_this()}, true);

    if (subsample > 1)
      pMapper = pMapper->superOrSubSampleDM(subsample, {multiLevelContext.getGridForLevel(0).shared_from_this()}, false);

  }

  void reflag() {

    if (pInputMapper != nullptr) {
      pMapper->unflagAllParticles();
      pInputMapper->flagParticles(flaggedParticles);
      pInputMapper->extendParticleListToUnreferencedGrids(multiLevelContext);
    }
    else {
      pMapper->unflagAllParticles();
      pMapper->flagParticles(flaggedParticles);
    }
  }


  virtual void write() {
    using namespace io;

    if(!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    initialiseParticleGenerator();

    cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
    cerr << (*pMapper);

    T boxlen = multiLevelContext.getGridForLevel(0).simsize;

    switch(outputFormat) {
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

    haveInitialisedRandomComponent=true;
  }



protected:

  int deepestLevelWithParticlesSelected() {
    for(size_t i=multiLevelContext.getNumLevels()-1; i>=0; --i) {
      if(multiLevelContext.getGridForLevel(i).hasFlaggedCells())
        return i;
    }
    throw std::runtime_error("No level has any particles selected");
  }

  int deepestLevel() {
    //TODO: can this be removed?
    return multiLevelContext.getNumLevels();
  }

  T get_wrapped_delta(T x0, T x1) {
    return multiLevelContext.getGridForLevel(0).getWrappedDelta(x0, x1);
  }


  void getCentre() {
    x0 = 0;
    y0 = 0;
    z0 = 0;

    int level = deepestLevelWithParticlesSelected();

    std::vector<size_t> particleArray;
    Grid<T> &grid = multiLevelContext.getGridForLevel(level);
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

  void centreParticle(long id) {
    std::tie(x0, y0, z0) = multiLevelContext.getGridForLevel(0).getCellCentroid(id);
  }

  void selectNearest() {
    auto & grid = multiLevelContext.getGridForLevel(deepestLevel()-1);
    pMapper->unflagAllParticles();
    size_t id = grid.getClosestIdNoWrap(Coordinate<T>(x0, y0, z0));
    cerr << "selectNearest " <<x0 << " " << y0 << " " << z0 << " " << id << " " << endl;
    grid.flagCells({id});

  }

  void select(std::function<bool(T, T, T)> inclusionFunction) {
    T delta_x, delta_y, delta_z;
    T xp, yp, zp;

    flaggedParticles.clear();



    for_each_level(level) {
      std::vector<size_t> particleArray;
      Grid<T> &grid = multiLevelContext.getGridForLevel(level);
      size_t N = grid.size3;
      for (size_t i = 0; i < N; i++) {
        std::tie(xp, yp, zp) = grid.getCellCentroid(i);
        delta_x = get_wrapped_delta(xp, x0);
        delta_y = get_wrapped_delta(yp, y0);
        delta_z = get_wrapped_delta(zp, z0);
        if (inclusionFunction(delta_x, delta_y, delta_z))
          particleArray.push_back(i);
      }

      grid.unflagAllCells();
      grid.flagCells(particleArray);

    }
  }

  void selectSphere(float radius) {
    T r2 = radius * radius;
    select([r2](T delta_x, T delta_y, T delta_z) -> bool {
      T r2_i = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
      return r2_i < r2;
    });

  }

  void selectCube(float side) {
    T side_by_2 = side / 2;
    select([side_by_2](T delta_x, T delta_y, T delta_z) -> bool {
      return abs(delta_x) < side_by_2 && abs(delta_y) < side_by_2 && abs(delta_z) < side_by_2;
    });
  }


  void setCentre(T xin, T yin, T zin) {
    x0 = xin;
    y0 = yin;
    z0 = zin;
  }

  auto calcConstraint(string name_in) {
    auto constraint = constraintGenerator.calcConstraintForAllLevels(name_in);
    constraint.toFourier();
    /*
    T norm = constraint.innerProduct(constraint).real();
    constraint/=sqrt(norm);
    constraint.convertToVector();

    auto constraint2 = constraintGenerator.calcConstraintForAllLevels(name_in);
    constraint2/=sqrt(norm);
    constraint2.toFourier();


    // constraint.setStateRecombined();

    constraint.toFourier();

    //constraint.getFieldForLevel(0).applyFilter(constraint.getFilterForLevel(0));
    //constraint.getFieldForLevel(1).applyFilter(constraint.getFilterForLevel(1));
    constraint.getFieldForLevel(2).applyFilter(constraint.getFilterForLevel(2));

    constraint.toReal();

    constraint.getFieldForLevel(2).addFieldFromDifferentGridWithFilter(constraint.getFieldForLevel(0),constraint.getFilterForLevel(0));
    constraint.getFieldForLevel(2).addFieldFromDifferentGrid(constraint.getFieldForLevel(1));


    constraint.toFourier();

    constraint.template setupFilters<MultiLevelRecombinedFilterFamily<T>>();


    cerr << "THING THING HERE HERE " << constraint2.innerProduct(constraint) << endl;
    constraint.toReal();
    for(size_t level=0; level<3; level++) {
      cerr << level << ": ";
      constraint.getFilterForLevel(level).debugInfo(cerr);
      cerr << endl;
      dumpGridData(level, constraint.getFieldForLevel(level).getDataVector());
    }

    assert(false);
    */
    constraint.toFourier();
    return constraint;
  }

  void calculate(string name) {
    if(!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    auto constraint_field = calcConstraint(name);
    auto val = constraint_field.innerProduct(outputField);

    cout << name << ": calculated value = " << val << endl;
  }

  virtual void constrain(string name, string type, float value) {
    bool relative = false;
    if (strcasecmp(type.c_str(), "relative") == 0) {
      relative = true;
    } else if (strcasecmp(type.c_str(), "absolute") != 0) {
      throw runtime_error("Constraint type must be either 'relative' or 'absolute'");
    }

    std::complex<T> constraint = value;
    auto vec = calcConstraint(name);

    std::complex<T> initv = vec.innerProduct(outputField);

    if (relative) constraint *= initv;

    cout << name << ": initial value = " << initv << ", constraining to " << constraint << endl;
    constraintApplicator.add_constraint(std::move(vec), constraint, initv);

  }

  void cov() {
    constraintApplicator.print_covariance();
  }


  virtual void fixConstraints() {
    if(!haveInitialisedRandomComponent)
      initialiseRandomComponent();

    constraintApplicator.prepare();
    constraintApplicator.applyConstraints();
  }

  virtual void done() {
    T pre_constraint_chi2 = outputField.getChi2();
    cerr << "BEFORE constraints chi^2=" << pre_constraint_chi2 << endl;
    fixConstraints();
    T post_constraint_chi2 = outputField.getChi2();
    cerr << "AFTER  constraints chi^2=" << post_constraint_chi2 << endl;
    cerr << "             delta-chi^2=" << post_constraint_chi2-pre_constraint_chi2 << endl;
    write();
  }

  void reverse() {
    for_each_level(level) {
      auto &field = outputField.getFieldForLevel(level);
      size_t N = field.getGrid().size3;
      for (size_t i = 0; i < N; i++)
        field[i] = -field[i];
    }
  }

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
      auto &fieldOriginal = fieldCopies[level];
      auto &field = outputField.getFieldForLevel(level);
      auto &grid = field.getGrid();
      field.toFourier();
      T k2;
      size_t N = grid.size3;
      for (size_t i = 0; i < N; i++) {
        k2 = grid.getFourierCellKSquared(i);
        if (k2 > k2max && k2 != 0) {
          field[i] = fieldOriginal[i];
        }
      }
    }

  }

  void reverseSmallK(T kmax) {

    T k2max = kmax * kmax;


    for_each_level(level) {
      T k2_g_min = std::numeric_limits<T>::max();
      T k2_g_max = 0.0;
      size_t modes_reversed = 0;
      auto &field = outputField.getFieldForLevel(level);
      field.toFourier();
      const Grid<T> &grid = field.getGrid();
      size_t tot_modes = grid.size3;
      T k2;
      size_t N = grid.size3;
      for (size_t i = 0; i < N; i++) {
        k2 = grid.getFourierCellKSquared(i);
        if (k2 < k2max && k2 != 0) {
          field[i] = -field[i];
          modes_reversed++;
        }
        if (k2 < k2_g_min && k2 != 0)
          k2_g_min = k2;
        if (k2 > k2_g_max)
          k2_g_max = k2;
      }
      cerr << "reverseSmallK: k reversal at " << sqrt(k2max) << "; grid was in range " << sqrt(k2_g_min) << " to " <<
      sqrt(k2_g_max) << endl;
      cerr << "               modes reversed = " << modes_reversed << " of " << tot_modes << endl;
    }

  }


};

#endif
