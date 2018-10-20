#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>


#include "tools/data_types/float_types.hpp"
#include "tools/numerics/vectormath.hpp"
#include "tools/numerics/fourier.hpp"
#include "tools/parser.hpp"
#include "io.hpp"
#include "dummyic.hpp"


#ifdef DOUBLEPRECISION
typedef double FloatType;
#else
typedef float FloatType;
#endif

// #define USE_COMPLEX_FFT
#ifdef USE_COMPLEX_FFT
using ICf = ICGenerator<complex<FloatType>> ;
#else
using ICf = ICGenerator<FloatType>;
#endif

void setup_parser(tools::ClassDispatch<ICf, void> &dispatch) {

  // Link command names used in parameter file to methods

  //Set the cosmology
  dispatch.add_class_route("TCMB", &ICf::setTCMB);
  dispatch.add_class_route("Om", &ICf::setOmegaM0);
  dispatch.add_class_route("Ob", &ICf::setOmegaB0);
  dispatch.add_class_route("Ol", &ICf::setOmegaLambda0);
  dispatch.add_class_route("hubble", &ICf::setHubble);
  dispatch.add_class_route("s8", &ICf::setSigma8);
  dispatch.add_class_route("ns", &ICf::setns);
  dispatch.add_class_route("zin", &ICf::setZ0);

  // Set seeds for random draws
  dispatch.add_class_route("seed", &ICf::setSeed);
  dispatch.add_class_route("seedfourier", &ICf::setSeedFourier);
  dispatch.add_class_route("seedfourier_parallel", &ICf::setSeedFourierParallel);
  dispatch.add_class_route("seedfourier_reverse", &ICf::setSeedFourierReverseOrder);

  // Optional computational properties
  dispatch.add_class_route("exact_power_spectrum_enforcement", &ICf::setExactPowerSpectrumEnforcement);
  dispatch.add_class_route("strays_on", &ICf::setStraysOn);
  dispatch.add_class_route("supersample", &ICf::setSupersample);
  dispatch.add_class_route("subsample", &ICf::setSubsample);

  // Grafic options
  dispatch.add_class_route("pvar", &ICf::setpvarValue);
  dispatch.add_class_route("center_grafic_output", &ICf::setCenteringOnRegion);

  // Gadget options
  dispatch.add_class_route("gadget_particle_type", &ICf::setGadgetParticleType);
  dispatch.add_class_route("gadget_flagged_particle_type", &ICf::setFlaggedGadgetParticleType);

  // Define input files
  dispatch.add_class_route("mapper_relative_to", &ICf::setInputMapper);
  dispatch.add_class_route("camb", &ICf::setCambDat);

  // Set output paths and format
  dispatch.add_class_route("outdir", &ICf::setOutDir);
  dispatch.add_class_route("outname", &ICf::setOutName);
  dispatch.add_class_route("outformat", &ICf::setOutputFormat);

  // Define grid structure
  dispatch.add_class_route("basegrid", &ICf::initBaseGrid);
  dispatch.add_class_route("zoomgrid", &ICf::initZoomGrid);
  dispatch.add_class_route("zoomgrid_with_origin_at", &ICf::initZoomGridWithLowLeftCornerAt);

  //Input Output of flagged particles
  dispatch.add_class_route("IDfile", &ICf::loadID);
  dispatch.add_class_route("append_IDfile", &ICf::appendID);
  dispatch.add_class_route("dump_IDfile", &ICf::dumpID);

  // Point to specific coordinates
  dispatch.add_class_route("centre_on", &ICf::centreParticle);
  dispatch.add_class_route("centre", &ICf::setCentre);

  // Flag cells
  dispatch.add_class_route("select_sphere", &ICf::selectSphere);
  dispatch.add_class_route("select_cube", &ICf::selectCube);
  dispatch.add_class_route("select_nearest", &ICf::selectNearest);
  dispatch.add_class_route("expand_flagged_region", &ICf::expandFlaggedRegion);
  dispatch.add_class_route("adapt_mask", &ICf::adaptMask);

  // Deal with modifications
  dispatch.add_class_route("calculate", &ICf::calculate);
  dispatch.add_class_route("filtering_scale", &ICf::setVarianceFilteringScale);
  dispatch.add_class_route("modify", &ICf::modify);
  dispatch.add_class_route("clear_modifications", &ICf::clearModifications);
  dispatch.add_class_route("done", &ICf::done);
  dispatch.add_class_route("apply_modifications", &ICf::applyModifications);
  dispatch.add_class_route("chi2", &ICf::getFieldChi2);

  dispatch.add_class_route("reverse", &ICf::reverse);
  dispatch.add_class_route("reverse_small_k", &ICf::reverseSmallK);
  dispatch.add_class_route("reseed_small_k", &ICf::reseedSmallK);

  // Write objects to files
  dispatch.add_class_route("dump_grid", &ICf::dumpGrid);
  dispatch.add_class_route("dump_ps", &ICf::dumpPS);
  dispatch.add_class_route("dump_tipsy", &ICf::saveTipsyArray);
  dispatch.add_class_route("dump_mask", &ICf::dumpMask);

  // Load existing random field instead of generating
  dispatch.add_class_route("import_level", &ICf::importLevel );

  // To debug, allow a level to be replaced with zeroes
  dispatch.add_class_route("zeroLevel", &ICf::zeroLevel);

  // Set an overall velocity offset for the entire box (useful for testing AMR sensitivity to flows)
  dispatch.add_class_route("velocity_offset", &ICf::setVelocityOffset);

}

void header(ostream &outf) {
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  outf << "# GM ICs code, compiled " << __DATE__ << " " << __TIME__ << endl;
  outf << "# git HEAD:" << GIT_VERSION << endl;
  if (sizeof(GIT_MODIFIED) > 1) {
    outf << "# However, the following files are modified:" << endl;
    outf << "#  " << GIT_MODIFIED << endl;
  }
  outf << "# Runtime: " << buf << endl << endl;
}

int main(int argc, char *argv[]) {
  using namespace std;

  if (argc != 2) {
    cerr << "Usage: ./IC paramfile" << endl;
    return -1;
  }

  ifstream inf;
  string fname(argv[1]);

  inf.open(fname);

  if (!inf.is_open()) {
    cerr << "Error: could not open parameter file " << argv[1] << endl;
    exit(1);
  }

  tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(fname));

  ofstream outf;
  outf.open("IC_output.params");

  if (!outf.is_open()) {
    cerr << "Error: could not open output file" << endl;
    exit(1);
  }

  header(outf);
  header(cerr);

  // Set up the command interpreter to issue commands to main_generator
  tools::ClassDispatch<ICf, void> dispatch_generator;
  setup_parser(dispatch_generator);

  // The main program is contained in this class:
  ICf generator(dispatch_generator);


  auto dispatch = dispatch_generator.specify_instance(generator);

  // Initialising FFTW is not strictly required (it will be done later if missed here) but, since it might print
  // some output, it's neater to do it straight away:
  tools::numerics::fourier::initialise();

  // Process commands
  dispatch.run_loop(inf, outf);

  return 0;
}
