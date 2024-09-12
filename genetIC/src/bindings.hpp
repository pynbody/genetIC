#ifndef _BINDINGS_HPP_INCLUDED
#define _BINDINGS_HPP_INCLUDED

#include "tools/parser.hpp"
#include "ic.hpp"

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
    dispatch.add_class_route("splice", &ICType::splice);

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

template<typename T>
void runInterpreter(T& target, const std::string &fname, std::ofstream &outf) {
  runInterpreter<T>(target, fname, &outf);
}

template<typename T>
void runInterpreter(T& target, const std::string &fname, std::ofstream *outf = nullptr) {
  tools::ClassDispatch<T, void> dispatch_generator;
  setup_parser<T>(dispatch_generator);

  auto dispatch = dispatch_generator.specify_instance(target);

  std::ifstream inf;
  inf.open(fname);

  if (!inf.is_open())
    throw std::runtime_error("Cannot open IC paramfile "+fname);

  tools::ChangeCwdWhileInScope temporaryCwd(tools::getDirectoryName(fname));

  if(outf != nullptr) {
    dispatch.run_loop(inf, *outf);
  } else {
    dispatch.run_loop(inf);
  }

}



#endif