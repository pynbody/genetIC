#ifndef IC_SWIFT_HPP
#define IC_SWIFT_HPP

#include "gadgethdf.hpp"

namespace io {
  namespace swift {
    template<typename GridDataType, typename OutputFloatType>
    class SwiftHDFOutput : public gadgethdf::GadgetHDFOutput<GridDataType, OutputFloatType> {
    public:
      SwiftHDFOutput(double Boxlength,
                     particle::mapper::ParticleMapper<GridDataType> &mapper,
                     particle::SpeciesToGeneratorMap<GridDataType> &generators,
                     const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
                     int nFiles) :
        gadgethdf::GadgetHDFOutput<GridDataType, OutputFloatType>(Boxlength, mapper, generators, cosmology, nFiles) {}

      bool forceVariableMass() const override {
        return true;
      }

    };

    //! \brief Generates swift output
    /*!
    \param name - name of output file
    \param Boxlength - simulation size in Mpc/h
    \param mapper - particle mapper used to link particles to grid locations
    \param generators - particles generators for each particle species (vector)
    \param cosmology - cosmological parameters
    */
    template<typename OutputFloatType, typename GridDataType>
    void save(const std::string &name, double Boxlength,
              particle::mapper::ParticleMapper<GridDataType> &mapper,
              particle::SpeciesToGeneratorMap<GridDataType> &generators,
              const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
              int nFiles) {

      HighFive::SilenceHDF5 silence; // suppresses HDF5 warnings while in scope
      SwiftHDFOutput<GridDataType, OutputFloatType> output(Boxlength, mapper, generators, cosmology, nFiles);
      output(name);
      logging::entry() << "Swift output saved.";
      logging::entry() << " Note that Swift output follows the GadgetHDF format and unit conventions.  You should";
      logging::entry() << " ensure that cleanup_h_factors and cleanup_velocity_factors are both enabled within your";
      logging::entry() << " Swift parameter file.";

    }

  }
}

#endif


