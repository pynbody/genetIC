#ifndef IC_TIPSY_HPP
#define IC_TIPSY_HPP

#include <src/tools/memmap.hpp>
#include "src/io.hpp"
#include "src/simulation/particles/mapper/mapper.hpp"
#include "src/simulation/particles/species.hpp"

namespace io {

  namespace tipsy {
    /*! \namespace io::tipsy
      \brief Functions related to outputting particles in tipsy format.
    */


    //! \struct io_header_tipsy
    /*! \brief Struct to hold tipsy header information.
    */
    struct io_header_tipsy {
      double scalefactor;
      int n;
      int ndim;
      int ngas;
      int ndark;
      int nstar;
    } header_tipsy;


    //! \brief Saves the specified multi-level field as a tipsy array.
    /*!
    \param filename - name of tipsy output file
    \param mapper - Mapper used to link particles to grid locations.
    \param generator - particle generator for this multi-level field.
    \param field - multi-level field to be output to tipsy format.
    */
    template<typename GridType, typename FloatType=tools::datatypes::strip_complex<GridType>>
    void saveFieldTipsyArray(const std::string &filename,
                             particle::mapper::ParticleMapper<GridType> &mapper,
                             particle::AbstractMultiLevelParticleGenerator<GridType> &generator,
                             fields::MultiLevelField<GridType> &field) {
      std::ofstream outfile(filename.c_str(), std::ofstream::binary);
      size_t lengthField = mapper.size();
      outfile.write(reinterpret_cast<char *>(&lengthField), 4);

      field.toReal();

      for (auto i = mapper.begin(generator); i != mapper.end(generator); ++i) {
        float data = float(tools::datatypes::real_part_if_complex(i.getField(field)));
        outfile.write(reinterpret_cast<char *>(&data), 4);
      }
    }

    namespace TipsyParticle {
      /*! \namespace io::tipsy::TipsyParticle
      \brief Classes related to defining tipsy particles.
      */

      //! \struct dark
      /*! \brief Dark matter particles.
      */
      struct dark {
        float mass, x, y, z, vx, vy, vz, eps, phi;
      };

      //! \struct gas
      /*! \brief Baryonic particles.
      */
      struct gas {
        float mass, x, y, z, vx, vy, vz, rho, temp, eps, metals, phi;
      };

      //! \brief Initialise dark matter particles
      template<typename T>
      void initialise(dark &p, const cosmology::CosmologicalParameters<T> & /*&cosmo*/) {
        p.phi = 0.0;
      }

      //! \brief Initialise baryonic particles:
      template<typename T>
      void initialise(gas &p, const cosmology::CosmologicalParameters<T> &cosmo) {
        p.temp = cosmo.TCMB * cosmo.scalefactorAtDecoupling / (cosmo.scalefactor * cosmo.scalefactor);
        p.metals = 0.0;
        p.rho = 0.0;
      }
    }


    //! \class TipsyOutput
    /*! \brief Class to handle outputting to tipsy format.
    */
    template<typename GridDataType, typename FloatType=tools::datatypes::strip_complex<GridDataType>>
    class TipsyOutput {
    protected:
      particle::SpeciesToGeneratorMap<GridDataType> generators; //!< Particle generators for each species.
      tools::MemMapFileWriter writer; //!< Writer used to process output file using memory maps.
      std::ofstream photogenic_file; //!< Photogenic output file.
      size_t iord; //!< Cumulative index offset.
      double pos_factor; //!< Factor to multiply internal position offset by to get tipsy units.
      double vel_factor; //!< Factor to multiply velocity offset by (especially as internal gadget-units velocities are used).
      double mass_factor;//!< Factor to multiply masses by to get tipsy units.
      double min_mass; //!< Minimum mass of particles in the simulation.
      double max_mass; //!< Maximum mass of particles in the simulation.
      double boxLength; //!< Simulation size in Mpc/h.
      std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pMapper; //!< Particle mapper.
      const cosmology::CosmologicalParameters<FloatType> &cosmology; //!< Cosmological paramters.

      //! \brief Save a block of tipsy particles in parallel
      template<typename ParticleType>
      void saveNextBlockOfTipsyParticles(particle::mapper::MapperIterator<GridDataType> &begin,
                                         size_t nMax = 256 * 1024 * 1024) {

        size_t n = std::min({nMax, begin.getNumRemainingParticles()});

        auto p = writer.getMemMap<ParticleType>(n);


        begin.parallelIterate([&](size_t i, const particle::mapper::MapperIterator<GridDataType> &localIterator) {
          TipsyParticle::initialise(p[i], cosmology);
          auto thisParticle = localIterator.getParticle();
          p[i].x = thisParticle.pos.x * pos_factor - 0.5;
          p[i].y = thisParticle.pos.y * pos_factor - 0.5;
          p[i].z = thisParticle.pos.z * pos_factor - 0.5;
          p[i].eps = thisParticle.soft * pos_factor;

          p[i].vx = thisParticle.vel.x * vel_factor;
          p[i].vy = thisParticle.vel.y * vel_factor;
          p[i].vz = thisParticle.vel.z * vel_factor;
          p[i].mass = thisParticle.mass * mass_factor;

#ifdef _OPENMP
          if (thisParticle.mass == min_mass && omp_get_thread_num() == 0)
#else
            if(thisParticle.mass==min_mass)
#endif
            photogenic_file << iord + i << std::endl;

        }, n);

        iord += n;


      }


      //! \brief Saves particles in tipsy format
      template<typename ParticleType>
      void saveTipsyParticles(particle::mapper::MapperIterator<GridDataType> &&begin,
                              particle::mapper::MapperIterator<GridDataType> &&end) {

        ParticleType p;
        TipsyParticle::initialise(p, cosmology);
        auto i = begin;
        std::thread writerThread;
        while (i != end) {
          saveNextBlockOfTipsyParticles<ParticleType>(i);
        }
      }

      //! \brief Save tipsy particles in a single thread
      template<typename ParticleType>
      void saveTipsyParticlesSingleThread(particle::mapper::MapperIterator<GridDataType> &&begin,
                                          particle::mapper::MapperIterator<GridDataType> &&end) {

        ParticleType p;
        FloatType x, y, z, vx, vy, vz, mass, eps;

        for (auto i = begin; i != end; ++i) {
          auto thisParticle = i.getParticle();

          p.x = thisParticle.pos.x * pos_factor - 0.5;
          p.y = thisParticle.pos.y * pos_factor - 0.5;
          p.z = thisParticle.pos.z * pos_factor - 0.5;
          p.eps = thisParticle.soft * pos_factor;

          p.vx = thisParticle.vel.x * vel_factor;
          p.vy = thisParticle.vel.y * vel_factor;
          p.vz = thisParticle.vel.z * vel_factor;
          p.mass = thisParticle.mass * mass_factor;

          writer.write<ParticleType>(p);

          if (mass == min_mass) {
            photogenic_file << iord << std::endl;
          }

          ++iord;
        }
      }


    public:
      //! \brief Constructor
      TipsyOutput(double boxLength,
                  const particle::SpeciesToGeneratorMap<GridDataType> &_generators,
                  std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pMapper,
                  const cosmology::CosmologicalParameters<FloatType> &cosmology) : generators(_generators),
                                                                                   iord(0),
                                                                                   boxLength(boxLength),
                                                                                   pMapper(pMapper),
                                                                                   cosmology(cosmology) {

      }

      //! \brief Main save operation
      void operator()(const std::string &filename) {

        // originally:
        // pmass in 1e10 h^-1 Msol
        // pos in Mpc h^-1
        // vel in km s^-1 a^1/2


        min_mass = std::numeric_limits<double>::max();
        max_mass = 0.0;

        FloatType mass, tot_mass = 0.0;

        const auto end = pMapper->end(*generators[particle::species::dm]); // don't want to keep re-evaluating this
        for (auto i = pMapper->begin(*generators[particle::species::dm]); i != end; ++i) {
          mass = i.getMass(); // sometimes can be MUCH faster than getParticle
          if (min_mass > mass) min_mass = mass;
          if (max_mass < mass) max_mass = mass;
          tot_mass += mass;
        }

        if (min_mass != max_mass) {
          photogenic_file.open("photogenic.txt");
        }

        mass_factor = cosmology.OmegaM0 / tot_mass; // tipsy convention: sum(mass)=Om0
        pos_factor = 1. / boxLength;              // boxsize = 1

        double dKpcUnit = boxLength * 1000 / cosmology.hubble;
        double dMsolUnit = 1e10 / cosmology.hubble / mass_factor;
        double dKmsUnit = sqrt(4.3022682e-6 * dMsolUnit / (dKpcUnit));

        vel_factor = std::pow(cosmology.scalefactor, -0.5) / dKmsUnit;


        io_header_tipsy header;

        header.scalefactor = cosmology.scalefactor;
        header.n = pMapper->size();
        header.ndim = 3;
        header.ngas = pMapper->size_gas();
        header.ndark = pMapper->size_dm();
        header.nstar = 0;


        std::ofstream paramfile;
        paramfile.open("tipsy.param");

        paramfile << "dKpcUnit = " << dKpcUnit << std::endl;
        paramfile << "dMsolUnit = " << dMsolUnit << std::endl;
        paramfile << "dHubble0 = " << 0.1 * cosmology.hubble * dKpcUnit / dKmsUnit << std::endl;
        paramfile << "bComove = 1 " << std::endl;

        paramfile.close();

        writer = tools::MemMapFileWriter(filename);

        writer.write<>(header);

        saveTipsyParticles<TipsyParticle::gas>(pMapper->beginGas(*generators.at(particle::species::baryon)),
                                               pMapper->endGas(*generators.at(particle::species::baryon)));
        saveTipsyParticles<TipsyParticle::dark>(pMapper->beginDm(*generators.at(particle::species::dm)),
                                                pMapper->endDm(*generators.at(particle::species::dm)));

      }
    };

    //! \brief Creates TipsyOutput object and uses it to save particles in tipsy format.
    /*!
    \param filename - name of tipsy output file
    \param boxLength - simulation size in Mpc/h
    \param generators - vector of particles generators for each particle type
    \param pMapper - Particle mapper, linking particles to grid locations.
    \param cosmology - cosmological parameters
    */
    template<typename GridDataType, typename T>
    void save(const std::string &filename, double boxLength,
              const particle::SpeciesToGeneratorMap<GridDataType> &generators,
              std::shared_ptr<particle::mapper::ParticleMapper<GridDataType>> pMapper,
              const cosmology::CosmologicalParameters<T> &cosmology) {

      TipsyOutput<GridDataType> output(boxLength, generators, pMapper, cosmology);
      output(filename);
    }


  }
}

#endif //IC_TIPSY_HPP
