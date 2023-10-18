#include <stdlib.h>
#include <sstream>
#include <ctime>
#include <iostream>
#include <iomanip>

#include "tools/parser.hpp"
#include "tools/logging.hpp"
#include "ic.hpp"
#include "dummyic.hpp"

int runGenetic(string paramFilename);

#ifdef USE_COMPLEX_FFT
using ICf = ICGenerator<complex<float>>;
using ICd = ICGenerator<complex<double>>;
#else
using ICf = ICGenerator<float>;
using ICd = ICGenerator<double>;
#endif

template<typename ICType>
void setup_parser(tools::ClassDispatch<ICType, void> &dispatch) {

  using FloatType = typename ICType::T;

  // Link command names used in parameter file to methods

  // Set the cosmology
  dispatch.add_class_route("TCMB", &ICType::setTCMB);
  dispatch.add_class_route("Om", &ICType::setOmegaM0);
  dispatch.add_class_route("Ob", &ICType::setOmegaB0);
  dispatch.add_class_route("Ol", &ICType::setOmegaLambda0);
  dispatch.add_class_route("hubble", &ICType::setHubble);
  dispatch.add_class_route("s8", &ICType::setSigma8);
  dispatch.add_class_route("ns", &ICType::setns);
  dispatch.add_class_route("zin", &ICType::setZ0);
  dispatch.add_class_route("z", &ICType::setZ0);

  // Set seeds for random draws
  // Static casts here needed to differentiate between overloaded versions of setSeed,
  // now that we have both DM and baryon fields to seed.
  // DEPRECATED forms:
  dispatch.add_deprecated_class_route("seed", "random_seed_real_space", static_cast<void (ICType::*)(int)>(&ICType::setSeed));
  dispatch.add_deprecated_class_route("seedfourier", "random_seed_serial", static_cast<void (ICType::*)(int)>(&ICType::setSeedFourier));
  dispatch.add_class_route("seedfourier_reverse", static_cast<void (ICType::*)(int)>(&ICType::setSeedFourierReverseOrder));
  dispatch.add_deprecated_class_route("seedfourier_parallel", "random_seed", static_cast<void (ICType::*)(int)>(&ICType::setSeedFourierParallel));

  // NEW FORMS, May 2020, ahead of release (to steer people towards a sensible default)
  dispatch.add_class_route("random_seed", static_cast<void (ICType::*)(int)>(&ICType::setSeedFourierParallel));
  dispatch.add_class_route("random_seed_serial", static_cast<void (ICType::*)(int)>(&ICType::setSeedFourier));
  dispatch.add_class_route("random_seed_real_space", static_cast<void (ICType::*)(int)>(&ICType::setSeed));

  // Optional computational properties
  dispatch.add_deprecated_class_route("exact_power_spectrum_enforcement", "fix_power", &ICType::setExactPowerSpectrumEnforcement);
  dispatch.add_class_route("fix_power", &ICType::setExactPowerSpectrumEnforcement);

  dispatch.add_class_route("strays_on", &ICType::setStraysOn);
  dispatch.add_class_route("supersample", &ICType::setSupersample);
  dispatch.add_class_route("supersample_gas", &ICType::setSupersampleGas);
  dispatch.add_class_route("subsample", &ICType::setSubsample);
  dispatch.add_class_route("eps_norm", &ICType::setEpsNorm);
  dispatch.add_class_route("num_neighbors", &ICType::setNumNeighbours);
  dispatch.add_class_route("num_neighbours", &ICType::setNumNeighbours);

  // Grafic options
  dispatch.add_class_route("pvar", &ICType::setpvarValue);
  dispatch.add_deprecated_class_route("center_grafic_output", "center_output",
                           &ICType::setCenteringOnRegion);


  dispatch.add_class_route("centre_output", &ICType::setCenteringOnRegion);
  // US variant:
  dispatch.add_class_route("center_output", &ICType::setCenteringOnRegion);

  // Gadget options
  dispatch.add_class_route("gadget_particle_type", &ICType::setGadgetParticleType);
  dispatch.add_class_route("gadget_flagged_particle_type", &ICType::setFlaggedGadgetParticleType);
  dispatch.add_class_route("gadget_num_files", &ICType::setGadgetNumFiles);

  // Define input files
  dispatch.add_class_route("mapper_relative_to", &ICType::setInputMapper);
  dispatch.add_class_route("camb", &ICType::setCambDat);
  dispatch.add_class_route("powerlaw_amplitude", &ICType::setPowerLawAmplitude);

  // Set output paths and format
  dispatch.add_class_route("outdir", &ICType::setOutDir);
  dispatch.add_class_route("outname", &ICType::setOutName);
  dispatch.add_class_route("outformat", &ICType::setOutputFormat);

  // Define grid structure - OLD NAMES
  dispatch.add_deprecated_class_route("basegrid", "base_grid", &ICType::initBaseGrid);
  dispatch.add_deprecated_class_route("zoomgrid", "zoom_grid", &ICType::initZoomGrid);
  dispatch.add_deprecated_class_route("zoomgrid_with_origin_at", "zoom_grid_with_origin_at", &ICType::initZoomGridWithLowLeftCornerAt);

  // Define grid structure - NEW NAMES May 2020
  dispatch.add_class_route("base_grid", &ICType::initBaseGrid);
  dispatch.add_class_route("zoom_grid", &ICType::initZoomGrid);
  dispatch.add_class_route("zoom_grid_with_origin_at", &ICType::initZoomGridWithLowLeftCornerAt);


  // Input Output of flagged particles - OLD NAMES
  dispatch.add_deprecated_class_route("IDfile", "id_file", &ICType::loadID);
  dispatch.add_deprecated_class_route("append_IDfile", "merge_id_file", &ICType::appendID);
  dispatch.add_deprecated_class_route("dump_IDfile","dump_id_file", &ICType::dumpID);

  // I/O of flagged particles - NEW NAMES May 2020
  dispatch.add_class_route("id_file", &ICType::loadID);
  dispatch.add_class_route("merge_id_file", &ICType::appendID);
  dispatch.add_class_route("dump_id_file", &ICType::dumpID);

  // Point to specific coordinates:
  dispatch.add_class_route("centre_on", &ICType::centreParticle);
  dispatch.add_class_route("centre", &ICType::setCentre);

  // US variants:
  dispatch.add_class_route("center_on", &ICType::centreParticle);
  dispatch.add_class_route("center", &ICType::setCentre);

  // Flag cells
  dispatch.add_class_route("select_sphere", &ICType::selectSphere);
  dispatch.add_class_route("select_ellipse", &ICType::selectEllipse);
  dispatch.add_class_route("select_cube", &ICType::selectCube);
  dispatch.add_class_route("select_nearest", &ICType::selectNearest);
  dispatch.add_class_route("expand_flagged_region", &ICType::expandFlaggedRegion);
  dispatch.add_class_route("adapt_mask", &ICType::adaptMask);
  dispatch.add_class_route("autopad", &ICType::setAutopad);

  // Deal with modifications
  dispatch.add_class_route("modify", &ICType::modify);
  dispatch.add_class_route("calculate", static_cast<void (ICType::*)(std::string)>(&ICType::calculate));
  dispatch.add_class_route("filtering_scale", &ICType::setVarianceFilteringScale);
  dispatch.add_class_route("clear_modifications", static_cast<void (ICType::*)()>(&ICType::clearModifications));
  dispatch.add_class_route("done", &ICType::done);
  dispatch.add_class_route("apply_modifications", static_cast<void (ICType::*)()>(&ICType::applyModifications));
  dispatch.add_class_route("chi2", static_cast<void (ICType::*)()>(&ICType::getFieldChi2));

  dispatch.add_class_route("reverse", static_cast<void (ICType::*)()>(&ICType::reverse));
  dispatch.add_class_route("reverse_small_k", static_cast<void (ICType::*)(FloatType)>(&ICType::reverseSmallK));
  dispatch.add_deprecated_class_route("splice_accuracy", "set_splice_accuracy", &ICType::set_splice_accuracy);
  dispatch.add_class_route("set_splice_accuracy", &ICType::set_splice_accuracy);
  dispatch.add_class_route("splice_seedfourier_parallel", &ICType::splice_seedfourier_parallel);
  dispatch.add_class_route("restart_splice", &ICType::restart_splice);
  dispatch.add_deprecated_class_route("splice", "splice_density", &ICType::splice_density);
  dispatch.add_class_route("splice_density", &ICType::splice_density);
  dispatch.add_class_route("splice_potential", &ICType::splice_potential);
  
  // Write objects to files
  // dispatch.add_class_route("dump_grid", &ICType::dumpGrid);
  dispatch.add_class_route("dump_grid", static_cast<void (ICType::*)(size_t)>(&ICType::dumpGrid));
  dispatch.add_class_route("dump_vx", static_cast<void (ICType::*)(size_t)>(&ICType::dumpVelocityX));
  dispatch.add_class_route("dump_grid_for_field", static_cast<void (ICType::*)(size_t, particle::species)>(&ICType::dumpGrid));
  dispatch.add_class_route("dump_grid_fourier", static_cast<void (ICType::*)(size_t)>(&ICType::dumpGridFourier));
  dispatch.add_class_route("dump_grid_fourier_for_field",
                           static_cast<void (ICType::*)(size_t, particle::species)>(&ICType::dumpGridFourier));
  dispatch.add_class_route("dump_ps", static_cast<void (ICType::*)(size_t)>(&ICType::dumpPS));
  dispatch.add_class_route("dump_ps_field", static_cast<void (ICType::*)(size_t, particle::species)>(&ICType::dumpPS));
  dispatch.add_class_route("dump_tipsy", static_cast<void (ICType::*)(std::string)>(&ICType::saveTipsyArray));
  dispatch.add_class_route("dump_tipsy_field", static_cast<void (ICType::*)(std::string, size_t)>(&ICType::saveTipsyArray));
  dispatch.add_class_route("dump_mask", &ICType::dumpMask);

  // Load existing random field instead of generating
  dispatch.add_class_route("import_level", static_cast<void (ICType::*)(size_t, std::string)>(&ICType::importLevel));

  // Extra commands related to the transfer functions:
  dispatch.add_class_route("baryon_tf_on", &ICType::setUsingBaryons);
  dispatch.add_class_route("baryons_all_levels", &ICType::setBaryonsOnAllLevels);

  // To debug
  dispatch.add_deprecated_class_route("zeroLevel", "zero_level", static_cast<void (ICType::*)(
    size_t)>(&ICType::zeroLevel)); // Retained for backwards compatibility.
  dispatch.add_class_route("zero_level", static_cast<void (ICType::*)(size_t)>(&ICType::zeroLevel)); // Use this instead.
  dispatch.add_class_route("zero_level_field", static_cast<void (ICType::*)(size_t, size_t)>(&ICType::zeroLevel));

  // Set an overall velocity offset for the entire box (useful for testing AMR sensitivity to flows)
  dispatch.add_class_route("velocity_offset", &ICType::setOffset);

}

void header(ostream &outf, std::string prefix="") {
  outf << prefix << "genetIC v1.5.0, compiled: " << __DATE__ << ", " << __TIME__ << endl;
  time_t now = time(0);
  struct tm tstruct;
  char buf[80];
  tstruct = *localtime(&now);
  strftime(buf, sizeof(buf), "%b %d %Y, %X", &tstruct);
  outf << prefix << "                 runtime: " << buf << endl << endl;
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


template<typename IC>
int runGenetic(string paramFilename) {
  ifstream inf;
  inf.open(paramFilename);

  if (!inf.is_open()) {
    logging::entry() << "Error: could not open parameter file " << paramFilename << endl;
    return -1;
  }

  tools::ChangeCwdWhileInScope temporary(tools::getDirectoryName(paramFilename));

  ofstream outf;
  outf.open("IC_output.params");

  if (!outf.is_open()) {
    logging::entry() << "Error: could not open output file" << endl;
    exit(1);
  }

  header(outf, "# ");
  header(cerr);

  // Set up the command interpreter to issue commands to main_generator
  tools::ClassDispatch<IC, void> dispatch_generator;
  setup_parser(dispatch_generator);

  // The main program is contained in this class:
  IC generator(dispatch_generator);


  auto dispatch = dispatch_generator.specify_instance(generator);

  // Initialising FFTW is not strictly required (it will be done later if missed here) but, since it might print
  // some output, it's neater to do it straight away:
  tools::numerics::fourier::initialise();

  // Process commands
  dispatch.run_loop(inf, outf);

  return 0;
}


void usageMessage() {
  using namespace std;
  cout << "Usage: genetIC paramfile [-f] [-c <cache-size>]" << endl << endl
       << " The paramfile is a text file of commands (see example provided with genetIC distribution)." << endl << endl
       << " If option -f is specified, genetIC uses float (32-bit) instead of double (64-bit) internally." << endl
       << " The output format is unaffected by the internal bit depth." << endl << endl
       << " If option -c <cache-size> is specified, the transfer function cache is limited to the specified"
          " number of entries. This can be used to reduce memory usage, but may slow down the code." << endl;
}

int main(int argc, char *argv[]) {
  using namespace std;

  bool useFloat = false;

  if (argc<2) {
    usageMessage();
    return -1;
  }


  string fname(argv[1]);

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "-f") == 0) {
      useFloat = true;
    } else if (strcmp(argv[i], "-c") == 0) {
      if (i + 1 >= argc) {
        cerr << "Error: -c option requires an argument" << endl;
        return -1;
      }
      cosmology::lru_cache_size = atoi(argv[++i]);
    } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      usageMessage();
      return 0;
    } else {
      if (i != argc - 1) {
        cerr << "Error: unknown option " << argv[i] << endl;
        return -1;
      }
    }
  }


  if(useFloat) {
    return runGenetic<ICf>(fname);
  } else {
    return runGenetic<ICd>(fname);
  }


}
