#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <complex>
#include <algorithm>
#include <iterator>
#include <sys/resource.h>

#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include "float_types.hpp"
#include "vectormath.hpp"
#include "fft.hpp"
#include "cosmo.hpp"
#include "parser.hpp"
#include "grid.hpp"
#include "multilevelcontext.hpp"
#include "src/particles/mapper.hpp"
#include "io.hpp"
#include "ic.hpp"
#include "dummyic.hpp"


#ifdef DOUBLEPRECISION
typedef double FloatType;
#else
typedef float FloatType;
#endif

#define USE_COMPLEX_FFT
#ifdef USE_COMPLEX_FFT
using ICf = ICGenerator<complex<FloatType>> ;
#else
using ICf = ICGenerator<FloatType>;
#endif

void setup_parser(ClassDispatch<ICf, void> &dispatch) {

  // Define the commands for the paramfile
  dispatch.add_class_route("TCMB", &ICf::setTCMB);
  dispatch.add_class_route("Om", &ICf::setOmegaM0);
  dispatch.add_class_route("Ob", &ICf::setOmegaB0);
  dispatch.add_class_route("Ol", &ICf::setOmegaLambda0);
  dispatch.add_class_route("hubble", &ICf::setHubble);
  dispatch.add_class_route("s8", &ICf::setSigma8);
  dispatch.add_class_route("basegrid", &ICf::initBaseGrid);
  dispatch.add_class_route("zoomgrid", &ICf::initZoomGrid);
  // dispatch.add_class_route("Boxl", &ICf::setBoxLen);
  dispatch.add_class_route("zin", &ICf::setZ0);
  // dispatch.add_class_route("n", &ICf::setn);
  dispatch.add_class_route("supersample", &ICf::setSupersample);
  dispatch.add_class_route("subsample", &ICf::setSubsample);
  dispatch.add_class_route("ns", &ICf::setns);
  dispatch.add_class_route("seed", &ICf::setSeed);
  dispatch.add_class_route("seedfourier", &ICf::setSeedFourier);
  dispatch.add_class_route("seedfourier_old", &ICf::setSeedFourierReverseOrder);
  dispatch.add_class_route("camb", &ICf::setCambDat);
  dispatch.add_class_route("outdir", &ICf::setOutDir);
  dispatch.add_class_route("outname", &ICf::setOutName);
  dispatch.add_class_route("outformat", &ICf::setOutputFormat);

  dispatch.add_class_route("IDfile", &ICf::loadID);
  dispatch.add_class_route("append_IDfile", &ICf::appendID);
  dispatch.add_class_route("dump_IDfile", &ICf::dumpID);
  dispatch.add_class_route("select_sphere", &ICf::selectSphere);
  dispatch.add_class_route("select_cube", &ICf::selectCube);
  dispatch.add_class_route("select_nearest", &ICf::selectNearest);
  dispatch.add_class_route("centre_on", &ICf::centreParticle);
  dispatch.add_class_route("centre", &ICf::setCentre);
  //dispatch.add_class_route("order", &ICf::reorderBuffer);
  // dispatch.add_class_route("truncate", &ICf::truncateBuffer);
  dispatch.add_class_route("calculate", &ICf::calculate);
  dispatch.add_class_route("constrain", &ICf::constrain);
  dispatch.add_class_route("cov", &ICf::cov);
  dispatch.add_class_route("done", &ICf::done);
  dispatch.add_class_route("fix_constraints", &ICf::fixConstraints);
  dispatch.add_class_route("offset_output", &ICf::offsetOutput);
  dispatch.add_class_route("reverse", &ICf::reverse);
  dispatch.add_class_route("reverse_small_k", &ICf::reverseSmallK);
  dispatch.add_class_route("reseed_small_k", &ICf::reseedSmallK);
  dispatch.add_class_route("dumpgrid", &ICf::dumpGrid);
  dispatch.add_class_route("dumpps", &ICf::dumpPS);
  dispatch.add_class_route("exact_power_spectrum_enforcement", &ICf::setExactPowerSpectrumEnforcement);

  dispatch.add_class_route("relative_to", &ICf::setInputMapper);
  // dispatch.add_class_route("zoom", &ICf::setZoom);
  // dispatch.add_class_route("n2", &ICf::setn2);
  // dispatch.add_class_route("zoom_IDfile", &ICf::setZoomParticles);
  // dispatch.add_class_route("dozoom", &ICf::doZoom);
  dispatch.add_class_route("dump", &ICf::saveTipsyArray);
  // dispatch.add_class_route("writeLevel", &ICf::writeLevel);

  dispatch.add_class_route("zeroLevel", &ICf::zeroLevel);

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

  ChangeCwdWhileInScope temporary(getDirectoryName(fname));

  ofstream outf;
  outf.open("IC_output.params");

  if (!outf.is_open()) {
    cerr << "Error: could not open output file" << endl;
    exit(1);
  }

  header(outf);
  header(cerr);

  // Set up the command interpreter to issue commands to main_generator
  ClassDispatch<ICf, void> dispatch_generator;
  setup_parser(dispatch_generator);

  // The main program is contained in this class:
  ICf generator(dispatch_generator);


  auto dispatch = dispatch_generator.specify_instance(generator);

  dispatch.run_loop(inf, outf);

  return 0;
}
