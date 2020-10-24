#include <stdlib.h>
#include <sstream>
#include <ctime>
#include <iostream>
#include <iomanip>

#include "tools/parser.hpp"
#include "tools/logging.hpp"
#include "ic.hpp"
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

  // Set the cosmology
  dispatch.add_class_route("TCMB", &ICf::setTCMB);
  dispatch.add_class_route("Om", &ICf::setOmegaM0);
  dispatch.add_class_route("Ob", &ICf::setOmegaB0);
  dispatch.add_class_route("Ol", &ICf::setOmegaLambda0);
  dispatch.add_class_route("hubble", &ICf::setHubble);
  dispatch.add_class_route("s8", &ICf::setSigma8);
  dispatch.add_class_route("ns", &ICf::setns);
  dispatch.add_class_route("zin", &ICf::setZ0);
  dispatch.add_class_route("z", &ICf::setZ0);

  // Set seeds for random draws
  // Static casts here needed to differentiate between overloaded versions of setSeed,
  // now that we have both DM and baryon fields to seed.
  // DEPRECATED forms:
  dispatch.add_deprecated_class_route("seed", "random_seed_real_space", static_cast<void (ICf::*)(int)>(&ICf::setSeed));
  dispatch.add_deprecated_class_route("seedfourier", "random_seed_serial", static_cast<void (ICf::*)(int)>(&ICf::setSeedFourier));
  dispatch.add_class_route("seedfourier_reverse", static_cast<void (ICf::*)(int)>(&ICf::setSeedFourierReverseOrder));
  dispatch.add_deprecated_class_route("seedfourier_parallel", "random_seed", static_cast<void (ICf::*)(int)>(&ICf::setSeedFourierParallel));

  // NEW FORMS, May 2020, ahead of release (to steer people towards a sensible default)
  dispatch.add_class_route("random_seed", static_cast<void (ICf::*)(int)>(&ICf::setSeedFourierParallel));
  dispatch.add_class_route("random_seed_serial", static_cast<void (ICf::*)(int)>(&ICf::setSeedFourier));
  dispatch.add_class_route("random_seed_real_space", static_cast<void (ICf::*)(int)>(&ICf::setSeed));

  // Optional computational properties
  dispatch.add_deprecated_class_route("exact_power_spectrum_enforcement", "fix_power", &ICf::setExactPowerSpectrumEnforcement);
  dispatch.add_class_route("fix_power", &ICf::setExactPowerSpectrumEnforcement);

  dispatch.add_class_route("strays_on", &ICf::setStraysOn);
  dispatch.add_class_route("supersample", &ICf::setSupersample);
  dispatch.add_class_route("supersample_gas", &ICf::setSupersampleGas);
  dispatch.add_class_route("subsample", &ICf::setSubsample);
  dispatch.add_class_route("eps_norm", &ICf::setEpsNorm);

  // Grafic options
  dispatch.add_class_route("pvar", &ICf::setpvarValue);
  dispatch.add_deprecated_class_route("center_grafic_output", "center_output",
                           &ICf::setCenteringOnRegion);


  dispatch.add_class_route("centre_output", &ICf::setCenteringOnRegion);
  // US variant:
  dispatch.add_class_route("center_output", &ICf::setCenteringOnRegion);

  // Gadget options
  dispatch.add_class_route("gadget_particle_type", &ICf::setGadgetParticleType);
  dispatch.add_class_route("gadget_flagged_particle_type", &ICf::setFlaggedGadgetParticleType);

  // Define input files
  dispatch.add_class_route("mapper_relative_to", &ICf::setInputMapper);
  dispatch.add_class_route("camb", &ICf::setCambDat);
  dispatch.add_class_route("powerlaw_amplitude", &ICf::setPowerLawAmplitude);

  // Set output paths and format
  dispatch.add_class_route("outdir", &ICf::setOutDir);
  dispatch.add_class_route("outname", &ICf::setOutName);
  dispatch.add_class_route("outformat", &ICf::setOutputFormat);

  // Define grid structure - OLD NAMES
  dispatch.add_deprecated_class_route("basegrid", "base_grid", &ICf::initBaseGrid);
  dispatch.add_deprecated_class_route("zoomgrid", "zoom_grid", &ICf::initZoomGrid);
  dispatch.add_deprecated_class_route("zoomgrid_with_origin_at", "zoom_grid_with_origin_at", &ICf::initZoomGridWithLowLeftCornerAt);

  // Define grid structure - NEW NAMES May 2020
  dispatch.add_class_route("base_grid", &ICf::initBaseGrid);
  dispatch.add_class_route("zoom_grid", &ICf::initZoomGrid);
  dispatch.add_class_route("zoom_grid_with_origin_at", &ICf::initZoomGridWithLowLeftCornerAt);


  // Input Output of flagged particles - OLD NAMES
  dispatch.add_deprecated_class_route("IDfile", "id_file", &ICf::loadID);
  dispatch.add_deprecated_class_route("append_IDfile", "merge_id_file", &ICf::appendID);
  dispatch.add_deprecated_class_route("dump_IDfile","dump_id_file", &ICf::dumpID);

  // I/O of flagged particles - NEW NAMES May 2020
  dispatch.add_class_route("id_file", &ICf::loadID);
  dispatch.add_class_route("merge_id_file", &ICf::appendID);
  dispatch.add_class_route("dump_id_file", &ICf::dumpID);

  // Point to specific coordinates:
  dispatch.add_class_route("centre_on", &ICf::centreParticle);
  dispatch.add_class_route("centre", &ICf::setCentre);

  // US variants:
  dispatch.add_class_route("center_on", &ICf::centreParticle);
  dispatch.add_class_route("center", &ICf::setCentre);

  // Flag cells
  dispatch.add_class_route("select_sphere", &ICf::selectSphere);
  dispatch.add_class_route("select_cube", &ICf::selectCube);
  dispatch.add_class_route("select_nearest", &ICf::selectNearest);
  dispatch.add_class_route("expand_flagged_region", &ICf::expandFlaggedRegion);
  dispatch.add_class_route("adapt_mask", &ICf::adaptMask);
  dispatch.add_class_route("autopad", &ICf::setAutopad);

  // Deal with modifications
  dispatch.add_class_route("modify", &ICf::modify);
  dispatch.add_class_route("calculate", static_cast<void (ICf::*)(std::string)>(&ICf::calculate));
  dispatch.add_class_route("filtering_scale", &ICf::setVarianceFilteringScale);
  dispatch.add_class_route("clear_modifications", static_cast<void (ICf::*)()>(&ICf::clearModifications));
  dispatch.add_class_route("done", &ICf::done);
  dispatch.add_class_route("apply_modifications", static_cast<void (ICf::*)()>(&ICf::applyModifications));
  dispatch.add_class_route("chi2", static_cast<void (ICf::*)()>(&ICf::getFieldChi2));

  dispatch.add_class_route("reverse", static_cast<void (ICf::*)()>(&ICf::reverse));
  dispatch.add_class_route("reverse_small_k", static_cast<void (ICf::*)(FloatType)>(&ICf::reverseSmallK));

  // Write objects to files
  // dispatch.add_class_route("dump_grid", &ICf::dumpGrid);
  dispatch.add_class_route("dump_grid", static_cast<void (ICf::*)(size_t)>(&ICf::dumpGrid));
  dispatch.add_class_route("dump_vx", static_cast<void (ICf::*)(size_t)>(&ICf::dumpVelocityX));
  dispatch.add_class_route("dump_grid_for_field", static_cast<void (ICf::*)(size_t, particle::species)>(&ICf::dumpGrid));
  dispatch.add_class_route("dump_grid_fourier", static_cast<void (ICf::*)(size_t)>(&ICf::dumpGridFourier));
  dispatch.add_class_route("dump_grid_fourier_for_field",
                           static_cast<void (ICf::*)(size_t, particle::species)>(&ICf::dumpGridFourier));
  dispatch.add_class_route("dump_ps", static_cast<void (ICf::*)(size_t)>(&ICf::dumpPS));
  dispatch.add_class_route("dump_ps_field", static_cast<void (ICf::*)(size_t, particle::species)>(&ICf::dumpPS));
  dispatch.add_class_route("dump_tipsy", static_cast<void (ICf::*)(std::string)>(&ICf::saveTipsyArray));
  dispatch.add_class_route("dump_tipsy_field", static_cast<void (ICf::*)(std::string, size_t)>(&ICf::saveTipsyArray));
  dispatch.add_class_route("dump_mask", &ICf::dumpMask);

  // Load existing random field instead of generating
  dispatch.add_class_route("import_level", static_cast<void (ICf::*)(size_t, std::string)>(&ICf::importLevel));

  // Extra commands related to the transfer functions:
  dispatch.add_class_route("baryon_tf_on", &ICf::setUsingBaryons);
  dispatch.add_class_route("baryons_all_levels", &ICf::setBaryonsOnAllLevels);

  // To debug
  dispatch.add_deprecated_class_route("zeroLevel", "zero_level", static_cast<void (ICf::*)(
    size_t)>(&ICf::zeroLevel)); // Retained for backwards compatibility.
  dispatch.add_class_route("zero_level", static_cast<void (ICf::*)(size_t)>(&ICf::zeroLevel)); // Use this instead.
  dispatch.add_class_route("zero_level_field", static_cast<void (ICf::*)(size_t, size_t)>(&ICf::zeroLevel));

  // Set an overall velocity offset for the entire box (useful for testing AMR sensitivity to flows)
  dispatch.add_class_route("velocity_offset", &ICf::setOffset);

}

void header(ostream &outf, std::string prefix="") {
  outf << prefix << "genetIC v1.1.1, compiled: " << __DATE__ << ", " << __TIME__ << endl;
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%b %d %Y, %X", &tstruct);
  outf << prefix << "              runtime:  " << buf << endl << endl;
#ifdef GIT_VERSION
  if (sizeof(GIT_VERSION) > 1)
    outf << prefix << "git HEAD:" << GIT_VERSION << endl;
  if (sizeof(GIT_MODIFIED) > 1) {
    outf << prefix << "At the time of compilation, the following files had been modified:" << endl;
    std::istringstream modifFiles(GIT_MODIFIED);
    while(!modifFiles.eof()) {
      std::string thisFile;
      modifFiles >> thisFile;
      outf << prefix << "  " << thisFile << endl;
    }
  }
#endif


}

int main(int argc, char *argv[]) {
  using namespace std;

  if (argc != 2) {
    logging::entry() << "usage: genetIC paramfile" << endl;
    return -1;
  }

  ifstream inf;
  string fname(argv[1]);

  inf.open(fname);

  if (!inf.is_open()) {
    logging::entry() << "Error: could not open parameter file " << argv[1] << endl;
    exit(1);
  }

  tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(fname));

  ofstream outf;
  outf.open("IC_output.params");

  if (!outf.is_open()) {
    logging::entry() << "Error: could not open output file" << endl;
    exit(1);
  }

  header(outf, "# ");
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
