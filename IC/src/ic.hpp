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

#define for_each_level(level) for(size_t level=0; level<2 && n[level]>0; level++)

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
  RandomFieldGenerator<T> randomFieldGenerator;

  CAMB<T> spectrum;

  // Everything about the grids:
  T boxlen[2], dx[2];      // the box length of each grid
  int n[2];                // number of grid divisions along one axis

  int supersample, subsample;               // DM supersampling to perform on zoom grid, and subsampling on base grid

  std::vector<GridPtrType> pGrid;       // the objects that help us relate points on the grid

  T x_off[2], y_off[2], z_off[2]; // x,y,z offsets for subgrids

  T xOffOutput, yOffOutput, zOffOutput;

  int zoomfac; // the zoom factor, i.e. the ratio dx/dx[1]



  int out;
  io::OutputFormat outputFormat;

  string indir, inname, base;


  bool prepared;
  bool exactPowerSpectrum;

  std::vector<size_t> genericParticleArray;
  std::vector<size_t> zoomParticleArray;


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
    n[0] = -1;
    n[1] = -1; // no subgrid by default
    boxlen[0] = -1;
    x_off[0] = y_off[0] = z_off[0] = 0;
    x_off[1] = y_off[1] = z_off[1] = 0;
    prepared = false;
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

  void setBoxLen(T in) {
    boxlen[0] = in;
    initGrid();
  }

  void setZ0(T in) {
    cosmology.redshift = in;
    cosmology.scalefactor = 1. / (cosmology.redshift + 1.);
  }

  void setn(int in) {
    n[0] = in;
    initGrid();
  }

  virtual void initGrid(unsigned int level = 0) {

    if (n[level] < 0 || boxlen[level] < 0)
      return;

    dx[level] = boxlen[level] / n[level];

    if (pGrid.size() != level)
      throw std::runtime_error("Trying to re-initialize a grid level");

    pGrid.push_back(std::make_shared<Grid<T>>(boxlen[0], n[level], dx[level],
                                                    x_off[level], y_off[level], z_off[level]));


    updateParticleMapper();
    updateMultilevelContext();

  }

  void setns(T in) {
    cosmology.ns = in;
  }

  void setn2(int in) {
    n[1] = in;
  }

  void setZoom(int in) {
    // open a subgrid which is the specified factor smaller than
    // the parent grid
    boxlen[1] = boxlen[0] / in;
    zoomfac = in;
  }

  void setZoomParticles(string fname) {
    appendParticleIdFile(fname);
    doZoom();
  }

  void doZoom() {
    if (n[1] == 0)
      throw (std::runtime_error("Set n2 before specifying the zoom particles"));

    if (boxlen[1] == 0)
      throw (std::runtime_error("Set the zoom factor before specifying the zoom particles"));

    pGrid[0]->gatherParticleList(zoomParticleArray);

    // find boundaries
    int x0, x1, y0, y1, z0, z1;
    int x, y, z;

    x0 = y0 = z0 = n[0];
    x1 = y1 = z1 = 0;
    Grid<T> g(n[0]);

    // TO DO: wrap the box sensibly

    for (size_t i = 0; i < zoomParticleArray.size(); i++) {
      std::tie(x, y, z) = g.getCellCoordinate(zoomParticleArray[i]);
      if (x < x0) x0 = x;
      if (y < y0) y0 = y;
      if (z < z0) z0 = z;
      if (x > x1) x1 = x;
      if (y > y1) y1 = y;
      if (z > z1) z1 = z;
    }

    // Now see if the zoom the user chose is OK
    int n_user = n[0] / zoomfac;
    if ((x1 - x0) > n_user || (y1 - y0) > n_user || (z1 - z0) > n_user) {
      throw (std::runtime_error(
        "Zoom particles do not fit in specified sub-box. Decrease zoom, or choose different particles. (NB wrapping not yet implemented)"));
    }
    // At this point we know things fit. All we need to do is choose
    // the correct offset to get the particles near the centre of the
    // zoom box.

    // Here is the bottom left of the box:
    x = (x0 + x1) / 2 - n[0] / (2 * zoomfac);
    y = (y0 + y1) / 2 - n[0] / (2 * zoomfac);
    z = (z0 + z1) / 2 - n[0] / (2 * zoomfac);

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (z < 0) z = 0;
    if (x > n[0] - n[0] / zoomfac) x = n[0] - n[0] / zoomfac;
    if (y > n[0] - n[0] / zoomfac) y = n[0] - n[0] / zoomfac;
    if (z > n[0] - n[0] / zoomfac) z = n[0] - n[0] / zoomfac;


    x_off[1] = x_off[0] + x * dx[0];
    y_off[1] = y_off[0] + y * dx[0];
    z_off[1] = z_off[0] + z * dx[0];


    initZoom();

  }

  void initZoom() {
    zoomfac = (boxlen[0] / n[0]) / (boxlen[1] / n[1]);
    dx[1] = boxlen[1] / n[1];
    cout << "Initialized a zoom region:" << endl;
    cout << "  Subbox length   = " << boxlen[1] << " Mpc/h" << endl;
    cout << "  n[1]            = " << n[1] << endl;
    cout << "  dx[1]           = " << dx[1] << endl;
    cout << "  Zoom factor     = " << zoomfac << endl;
    cout << "  Low-left corner = " << x_off[1] << ", " << y_off[1] << ", " << z_off[1] << endl;
    cout << "  Num particles   = " << zoomParticleArray.size() << endl;
    initGrid(1);

    updateParticleMapper();

    cout << "  Total particles = " << pMapper->size() << endl;
  }


  void setOutputMode(int in) {
    out = in; // the joy of writing this line is substantial

    if (out != 0 && out != 1 && out != 2)
      throw runtime_error("Wrong output format, choose 0 (HDF5), 1 (Gadget) or 2 (both)");

#ifndef HAVE_HDF5
    if (out != 1)
      throw runtime_error("Not compiled with HDF5. Only output=1 is allowed in this case!");
#endif

    updateParticleMapper();

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
    updateMultilevelContext();
  }

  void setOutDir(std::string in) {
    indir = in;
  }

  void setOutName(std::string in) {
    inname = in;
  }

  void setGadgetFormat(int in) {
    outputFormat = static_cast<io::OutputFormat>(in);
    updateParticleMapper();
    if (!prepared)
      prepare(); //compatibility with old paramfiles
  }

  string make_base(string basename, int level = 0) {
    ostringstream nult;
    if (inname.size() == 0) {
      nult << basename << "IC_iter_" << floatinfo<T>::name << "_z" << cosmology.redshift << "_" << n[level] <<
      "_L" << boxlen[level];

    } else {
      if (level == 0)
        nult << basename << "/" << inname;
      else
        nult << basename << "/" << inname << "_" << level;
    }
    return nult.str();
  }

  void readCamb() {

    updateMultilevelContext();
  }

  void updateMultilevelContext() {
    if (spectrum.isUsable()) {
      multiLevelContext.clear();
      for_each_level(level) {
        multiLevelContext.addLevel(spectrum.getPowerSpectrumForGrid(*(this->pGrid[level])),
                                   this->pGrid[level]);
      }
    }
  }


  virtual void drawRandom() {
    randomFieldGenerator.draw();
  }

  virtual void zeroLevel(int level) {
    cerr << "*** Warning: your script calls zeroLevel(" << level << "). This is intended for testing purposes only!" <<
    endl;
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
    assert(data.size() == pGrid[level]->size3);
    ostringstream filename;
    filename << indir << "/grid-" << level << ".npy";

    const int dim[3] = {n[level], n[level], n[level]};
    numpy::SaveArrayAsNumpy(filename.str(), false, 3, dim, data.data());

    filename.str("");
    filename << indir << "/grid-info-" << level << ".txt";

    ofstream ifile;
    ifile.open(filename.str());
    cerr << "Writing to " << filename.str() << endl;

    ifile << x_off[level] << " " << y_off[level] << " " << z_off[level] << " " << boxlen[level] << endl;
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
    powsp_noJing(n[level], field.getDataVector(),
                 multiLevelContext.getCovariance(level),
                 (base + "_" + ((char) (level + '0')) + ".ps").c_str(), boxlen[level]);
  }


  virtual void zeldovich() {
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
    auto gridForOutput = pGrid[level];

    if (xOffOutput != 0 || yOffOutput != 0 || zOffOutput != 0) {
      gridForOutput = std::make_shared<OffsetGrid<T>>(pGrid[0], xOffOutput, yOffOutput, zOffOutput);
    }
    return gridForOutput;
  }

  void updateParticleMapper() {

    if (pGrid.size() == 0)
      return;

    // make a basic mapper for the base level grid
    pMapper = std::shared_ptr<particle::ParticleMapper<T>>(new particle::OneLevelParticleMapper<T>(
      getGridWithOutputOffset(0)
    ));


    if (pGrid.size() >= 3) {
      // possible future enhancement, but for now...
      throw runtime_error("Don't know how to set up a mapper for more than one level of refinement");
    }

    if (pGrid.size() == 2) {
      // it's a zoom!
      auto pMapperLevel1 = std::shared_ptr<particle::ParticleMapper<T>>(
        new particle::OneLevelParticleMapper<T>(getGridWithOutputOffset(1)));

      pMapper = std::shared_ptr<particle::ParticleMapper<T>>(
        new particle::TwoLevelParticleMapper<T>(pMapper, pMapperLevel1, zoomParticleArray, zoomfac * zoomfac * zoomfac));
    }

    if (cosmology.OmegaBaryons0 > 0) {

      // Add gas only to the deepest level. Pass the whole pGrid
      // vector if you want to add gas to every level.
      auto gasMapper = pMapper->addGas(cosmology.OmegaBaryons0 / cosmology.OmegaM0,
                                       {pGrid.back()});

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
      pMapper = pMapper->superOrSubSampleDM(supersample, {pGrid.back()}, true);

    if (subsample > 1)
      pMapper = pMapper->superOrSubSampleDM(subsample, {pGrid[0]}, false);

  }

  void clearAndDistributeParticleList() {

    if (pInputMapper != nullptr) {
      pMapper->clearParticleList();
      pInputMapper->distributeParticleList(genericParticleArray);
      pInputMapper->extendParticleListToUnreferencedGrids(this->pGrid);
    }
    else {
      pMapper->clearParticleList();
      pMapper->distributeParticleList(genericParticleArray);
    }
  }


  virtual void write() {
    using namespace io;

    zeldovich();

    cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
    cerr << (*pMapper);

    switch(outputFormat) {
      case OutputFormat::gadget2:
      case OutputFormat::gadget3:
        gadget::save(base + ".gadget", boxlen[0], *pMapper,
                     *pParticleGenerator,
                     cosmology, static_cast<int>(outputFormat));
        break;
      case OutputFormat::tipsy:
        tipsy::save(base + ".tipsy", boxlen[0], *pParticleGenerator,
                    pMapper, cosmology);
        break;
      case OutputFormat::grafic:
        grafic::save(base+".grafic", *pParticleGenerator, multiLevelContext, cosmology);
        break;
      default:
        throw std::runtime_error("Unknown output format");
    }

  }

  void makeInitialRealizationWithoutConstraints() {
    
    base = make_base(indir);

    drawRandom();

    applyPowerSpec();
  }


  virtual void prepare() {
    if (prepared)
      throw (std::runtime_error("Called prepare, but grid is already prepared for constraints"));

    makeInitialRealizationWithoutConstraints();

    prepared = true;

  }


protected:

  int deepestLevelWithParticlesSelected() {
    if (pGrid.size() > 1 && pGrid[1]->estimateParticleListSize() > 0) return 1; else return 0;
  }

  int deepestLevel() {
    if (pGrid.size() > 1) return 1; else return 0;
  }

  T get_wrapped_delta(T x0, T x1) {
    return pGrid[0]->getWrappedDelta(x0, x1);
  }


  void getCentre() {
    x0 = 0;
    y0 = 0;
    z0 = 0;

    int level = deepestLevelWithParticlesSelected();

    std::vector<size_t> particleArray;
    pGrid[level]->gatherParticleList(particleArray);

    auto p0_location = pGrid[level]->getCellCentroid(particleArray[0]);

    for (size_t i = 0; i < particleArray.size(); i++) {
      auto pi_location = pGrid[level]->getCellCentroid(particleArray[i]);
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

    io::getBuffer(genericParticleArray, filename);
    size_t size = genericParticleArray.size();
    std::sort(genericParticleArray.begin(), genericParticleArray.end());
    genericParticleArray.erase(std::unique(genericParticleArray.begin(), genericParticleArray.end()),
                               genericParticleArray.end());
    if (genericParticleArray.size() < size)
      cerr << "  ... erased " << size - genericParticleArray.size() << " duplicate particles" << endl;
    cerr << "  -> total number of particles is " << genericParticleArray.size() << endl;

    clearAndDistributeParticleList();
  }

  void loadParticleIdFile(std::string filename) {
    genericParticleArray.clear();
    appendParticleIdFile(filename);
  }


  auto calcConstraint(string name_in) {
    auto constraint = constraintGenerator.calcConstraintForAllLevels(name_in);
    constraint.toFourier();
    return constraint;
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

  void dumpID(string fname) {
    std::vector<size_t> results;
    cerr << "dumpID using current mapper:" << endl;
    cerr << (*pMapper);
    pMapper->gatherParticleList(results);
    io::dumpBuffer(results, fname);
  }

  void centreParticle(long id) {
    std::tie(x0, y0, z0) = pGrid[0]->getCellCentroid(id);
  }

  void selectNearest() {
    auto grid = pGrid[deepestLevel()];
    pMapper->clearParticleList();
    size_t id = grid->getClosestIdNoWrap(Coordinate<T>(x0, y0, z0));
    cerr << "selectNearest " <<x0 << " " << y0 << " " << z0 << " " << id << " " << endl;
    grid->distributeParticleList({id});

  }

  void select(std::function<bool(T, T, T)> inclusionFunction) {
    T delta_x, delta_y, delta_z;
    T xp, yp, zp;

    genericParticleArray.clear();



    for_each_level(level) {
      std::vector<size_t> particleArray;
      size_t N = this->pGrid[level]->size3;
      for (size_t i = 0; i < N; i++) {
        std::tie(xp, yp, zp) = pGrid[level]->getCellCentroid(i);
        delta_x = get_wrapped_delta(xp, x0);
        delta_y = get_wrapped_delta(yp, y0);
        delta_z = get_wrapped_delta(zp, z0);
        if (inclusionFunction(delta_x, delta_y, delta_z))
          particleArray.push_back(i);
      }

      pGrid[level]->clearParticleList();
      pGrid[level]->distributeParticleList(particleArray);

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

  void reorderBuffer() {

    cout << "Reordering buffer radially..." << endl;
    cout << " [taking centre = " << x0 << " " << y0 << " " << z0 << "]" << endl;

    std::vector<T> r2(genericParticleArray.size());
    T delta_x, delta_y, delta_z;

    std::vector<size_t> index(genericParticleArray.size());

    for (size_t i = 0; i < genericParticleArray.size(); i++) {
      std::tie(delta_x, delta_y, delta_z) = pGrid[0]->getCellCentroid(i);
      delta_x = get_wrapped_delta(delta_x, x0);
      delta_y = get_wrapped_delta(delta_y, y0);
      delta_z = get_wrapped_delta(delta_z, z0);
      r2[i] = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
      index[i] = i;
    }

    // Now sort the index array
    std::sort(index.begin(), index.end(),
              [&r2](size_t i1, size_t i2) { return r2[i1] < r2[i2]; });

    // Turn the index array into something pointing to the particles
    for (size_t i = 0; i < genericParticleArray.size(); i++) {
      index[i] = genericParticleArray[index[i]];
    }

    // Copy back into the particle array
    for (size_t i = 0; i < genericParticleArray.size(); i++) {
      genericParticleArray[i] = index[i];
    }

  }

  void truncateBuffer(float x) {
    if (x < 0) {
      cerr << "Truncate command takes a fraction between 0 and 1 or a number>=1" << endl;
      exit(0);
    }
    if (x < 1)
      genericParticleArray.resize(((int) (genericParticleArray.size() * x)));
    else
      genericParticleArray.resize(((int) x));
  }

  void calculate(string name) {
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
    constraintApplicator.prepare();
    constraintApplicator.applyConstraints();
  }

  virtual void done() {
    T pre_constraint_chi2 = outputField.getChi2();
    cerr << "BEFORE constraints: chi^2=" << pre_constraint_chi2 << endl;
    fixConstraints();
    T post_constraint_chi2 = outputField.getChi2();
    cerr << "AFTER  constraints: chi^2=" << post_constraint_chi2 << endl;
    cerr << "              delta chi^2=" << post_constraint_chi2-pre_constraint_chi2 << endl;
    write();
  }

  void reverse() {
    for_each_level(level) {
      auto &field = outputField.getFieldForLevel(level);
      size_t N = pGrid[level]->size3;
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
    makeInitialRealizationWithoutConstraints();

    // copy back the old field
    for_each_level(level) {
      auto &fieldOriginal = fieldCopies[level];
      auto &field = outputField.getFieldForLevel(level);
      field.toFourier();
      const auto &grid = *(this->pGrid[level]);
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
      size_t tot_modes = pGrid[level]->size3;
      auto &field = outputField.getFieldForLevel(level);
      field.toFourier();
      const Grid<T> &grid = *(this->pGrid[level]);
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
