#ifndef IC_GADGET_HPP
#define IC_GADGET_HPP

#include "src/tools/memmap.hpp"
#include "src/tools/data_types/float_types.hpp"
#include "src/io.hpp"
#include "src/simulation/particles/species.hpp"
#include <vector>

namespace io {
  /*!
  \namespace io::gadget
  \brief Classes related to outputting particle data in the gadget 2 and gadget 3 formats.
*/
  namespace gadget {


    using std::cerr;
    using std::endl;
    using std::vector;


    /*! \struct io_header_2
        \brief The data structure for the header of gadget2 format files
    */
    struct io_header_2 // header for gadget2
    {
      int npart[6];
      double mass[6];
      double time;
      double redshift;
      int flag_sfr;
      int flag_feedback;
      unsigned int nPartTotal[6];
      int flag_cooling;
      int num_files;
      double BoxSize;
      double Omega0;
      double OmegaLambda;
      double HubbleParam;
      // Also looks like we need the following, but we don't actually need to set these:
      int flag_stellarage;
      int flag_metals;
      // If more than 10^32 particles, this part stores the number of particles divided by 2^32 (rounding down) and nPartTotal contains the remainder.
      unsigned int nPartTotalHighWord[6];
      int flag_entropy_instead_u;
      char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 6 * 4 - 3 * 4];  /* fills to 256 Bytes */
    };


    /*! \struct io_header_3
        \brief The data staructure for the header of gadget3 format files
    */
    struct io_header_3 //header for gadget3
    {
      int npart[6];
      /*!< number of particles of each type in this file */
      double mass[6];
      /*!< mass of particles of each type. If 0, then the masses are explicitly
                     stored in the mass-block of the snapshot file, otherwise they are omitted */
      double time;
      /*!< time of snapshot file */
      double redshift;
      /*!< redshift of snapshot file */
      int flag_sfr;
      /*!< flags whether the simulation was including star formation */
      int flag_feedback;
      /*!< flags whether feedback was included (obsolete) */
      unsigned int nPartTotal[6];
      /*!< total number of particles of each type in this snapshot. This can be
                       different from npart if one is dealing with a multi-file snapshot. */
      int flag_cooling;
      /*!< flags whether cooling was included  */
      int num_files;
      /*!< number of files in multi-file snapshot */
      double BoxSize;
      /*!< box-size of simulation in case periodic boundaries were used */
      double Omega0;
      /*!< matter density in units of critical density */
      double OmegaLambda;
      /*!< cosmological constant parameter */
      double HubbleParam;
      /*!< Hubble parameter in units of 100 km/sec/Mpc */
      int flag_stellarage;
      /*!< flags whether the file contains formation times of star particles */
      int flag_metals;
      /*!< flags whether the file contains metallicity values for gas and star
                     particles */
      unsigned int nPartTotalHighWord[6];
      /*!< High word of the total number of particles of each type (see header2)*/
      int flag_entropy_instead_u;
      /*!< flags that IC-file contains entropy instead of u */
      int flag_doubleprecision;
      /*!< flags that snapshot contains double-precision instead of single precision */

      int flag_ic_info;
      /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                              or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                              For snapshots files, the value informs whether the simulation was evolved from
                              Zeldoch or 2lpt ICs. Encoding is as follows:
                                 FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                 FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                 FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                 FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                 FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                              All other values, including 0 are interpreted as "don't know" for backwards compatability.
                          */
      float lpt_scalingfactor;
      /*!< scaling factor for 2lpt initial conditions */

      char fill[48];    /*!< fills to 256 Bytes */

    };


    template<typename OutputFloatType, typename InternalFloatType>
    io_header_2 createGadget2Header(vector<InternalFloatType> masses, vector<size_t> npart,
                                    vector<size_t> npartTotal, int nFiles, double Boxlength,
                                    const cosmology::CosmologicalParameters<InternalFloatType> &cosmology) {
      io_header_2 header2;
      ::memset(&header2, 0, sizeof(io_header_2)); // ensure unused flags are all zero
      header2.npart[0] = npart[0];
      header2.npart[1] = npart[1];
      header2.npart[2] = npart[2];
      header2.npart[3] = npart[3];
      header2.npart[4] = npart[4];
      header2.npart[5] = npart[5];
      header2.mass[0] = masses[0];
      header2.mass[1] = masses[1];
      header2.mass[2] = masses[2];
      header2.mass[3] = masses[3];
      header2.mass[4] = masses[4];
      header2.mass[5] = masses[5];
      header2.time = cosmology.scalefactor;
      header2.redshift = cosmology.redshift;
      header2.flag_sfr = 0;
      header2.flag_feedback = 0;
      header2.nPartTotal[0] = (unsigned int) (npartTotal[0]);
      header2.nPartTotal[1] = (unsigned int) (npartTotal[1]);
      header2.nPartTotal[2] = (unsigned int) (npartTotal[2]);
      header2.nPartTotal[3] = (unsigned int) (npartTotal[3]);
      header2.nPartTotal[4] = (unsigned int) (npartTotal[4]);
      header2.nPartTotal[5] = (unsigned int) (npartTotal[5]);
      // Same basic thing should happen here as with gadget3:
      header2.nPartTotalHighWord[0] = (unsigned int) (npartTotal[0] >> 32);
      header2.nPartTotalHighWord[1] = (unsigned int) (npartTotal[1] >> 32);
      header2.nPartTotalHighWord[2] = (unsigned int) (npartTotal[2] >> 32);
      header2.nPartTotalHighWord[3] = (unsigned int) (npartTotal[3] >> 32);
      header2.nPartTotalHighWord[4] = (unsigned int) (npartTotal[4] >> 32);
      header2.nPartTotalHighWord[5] = (unsigned int) (npartTotal[5] >> 32);
      header2.flag_cooling = 0;
      header2.num_files = nFiles;
      header2.BoxSize = Boxlength;
      header2.Omega0 = cosmology.OmegaM0;
      header2.OmegaLambda = cosmology.OmegaLambda0;
      header2.HubbleParam = cosmology.hubble;

      // I don't think we use any of the extra flags:
      header2.flag_stellarage = 0;
      header2.flag_metals = 0;
      header2.flag_entropy_instead_u = 0;

      if (npartTotal[0] > 0) { //options for baryons
        header2.flag_sfr = 1;
        header2.flag_feedback = 1;
        header2.flag_cooling = 1;
      }


      return header2;
    }

    template<typename OutputFloatType, typename InternalFloatType>
    io_header_3 createGadget3Header(vector<InternalFloatType> masses, vector<size_t> npart,
                                    vector<size_t> npartTotal, int nFiles,
                                    double Boxlength,
                                    const cosmology::CosmologicalParameters<InternalFloatType> &cosmology) {
      io_header_3 header3;
      ::memset(&header3, 0, sizeof(io_header_3)); // ensure unused flags are all zero
      header3.npart[0] = (unsigned int) (npart[0]);
      header3.npart[1] = (unsigned int) (npart[1]);
      header3.npart[2] = (unsigned int) (npart[2]);
      header3.npart[3] = (unsigned int) (npart[3]);
      header3.npart[4] = (unsigned int) (npart[4]);
      header3.npart[5] = (unsigned int) (npart[5]);
      header3.mass[0] = (double) (masses[0]);
      header3.mass[1] = (double) (masses[1]);
      header3.mass[2] = (double) (masses[2]);
      header3.mass[3] = (double) (masses[3]);
      header3.mass[4] = (double) (masses[4]);
      header3.mass[5] = (double) (masses[5]);
      header3.time = cosmology.scalefactor;
      header3.redshift = cosmology.redshift;
      header3.flag_sfr = 0;
      header3.flag_feedback = 0;
      header3.nPartTotal[0] = (unsigned int) (npartTotal[0]);
      header3.nPartTotal[1] = (unsigned int) (npartTotal[1]);
      header3.nPartTotal[2] = (unsigned int) (npartTotal[2]);
      header3.nPartTotal[3] = (unsigned int) (npartTotal[3]);
      header3.nPartTotal[4] = (unsigned int) (npartTotal[4]);
      header3.nPartTotal[5] = (unsigned int) (npartTotal[5]);
      header3.flag_cooling = 0;
      header3.num_files = nFiles;
      header3.BoxSize = Boxlength;
      header3.Omega0 = cosmology.OmegaM0;
      header3.OmegaLambda = cosmology.OmegaLambda0;
      header3.HubbleParam = cosmology.hubble;
      header3.flag_stellarage = 0;  /*!< flags whether the file contains formation times of star particles */
      header3.flag_metals = 0;    /*!< flags whether the file contains metallicity values for gas and star  particles */
      header3.nPartTotalHighWord[0] = (unsigned int) (npartTotal[0] >> 32);
      header3.nPartTotalHighWord[1] = (unsigned int) (npartTotal[1] >> 32); //copied from Gadget3
      header3.nPartTotalHighWord[2] = (unsigned int) (npartTotal[2] >> 32);
      header3.nPartTotalHighWord[3] = (unsigned int) (npartTotal[3] >> 32);
      header3.nPartTotalHighWord[4] = (unsigned int) (npartTotal[4] >> 32);
      header3.nPartTotalHighWord[5] = (unsigned int) (npartTotal[5] >> 32);
      header3.flag_entropy_instead_u = 0; /*!< flags that IC-file contains entropy instead of u */
      header3.flag_doubleprecision = tools::datatypes::floatinfo<OutputFloatType>::doubleprecision;
      header3.flag_ic_info = 1;
      header3.lpt_scalingfactor = 0.f; /*!dummy value since we never use ic_info!=1 */

      if (npartTotal[0] > 0) { //options for baryons & special behavior
        header3.flag_sfr = 1;
        header3.flag_feedback = 1;
        header3.flag_cooling = 1;
        //none of the following is currently supported in this code:
        //header3.flag_stellarage=1;  /*!< flags whether the file contains formation times of star particles */
        //header3.flag_metals=1;    /*!< flags whether the file contains metallicity values for gas and star  particles */
        //header3.flag_entropy_instead_u=1; /*!< flags that IC-file contains entropy instead of u */
      }

      return header3;
    }


    /*! \class GadgetOutput
    \brief Class to handle output to gadget files.
    */
    template<typename GridDataType, typename OutputFloatType>
    class GadgetOutput {
    protected:
      using InternalFloatType = tools::datatypes::strip_complex<GridDataType>;

      particle::mapper::ParticleMapper<GridDataType> &mapper; //!< Particle mapper, for relating offsets in the file to GenetIC grid cells.
      particle::SpeciesToGeneratorMap<GridDataType> generators; //!< Particle generators for each particle species.
      const cosmology::CosmologicalParameters<InternalFloatType> &cosmology; //!< Struct containing cosmological parameters.
      std::vector<tools::MemMapFileWriter> writers; //!< Low-level file operations are handled by this object.

      size_t nTotal; //!< Total number of particles to output across all files and types
      std::vector<size_t> nPartPerFile; //!< Number of particles to write per file
      vector<size_t> nPartPerType; //!< Number of particles of each gadget type
      std::vector<std::vector<size_t>> nPartPerTypePerFile; //!< nPartPerTypePerFile[i][j] = num particles in file i of type j

      double boxLength; //!< Size of simulation box.
      int gadgetVersion; //!< Which version of the gadget file to output. Allowed values 2 or 3.
      vector<InternalFloatType> masses; //!< Masses of particles if constant. Zero if variable.

      bool variableMass; //!< Stores whether we are using variable mass gadget particles
      unsigned int nFiles = 1; //!< Number of files to write

      // Mapping between gadget particle types (0->6) and our internal field type. This selects the appropriate
      // transfer function, if multiple are being used.
      std::vector<particle::species> gadgetTypeToSpecies{particle::species::baryon,
                                                           particle::species::dm,
                                                           particle::species::dm,
                                                           particle::species::dm,
                                                           particle::species::dm,
                                                           particle::species::dm};

      //! \brief Save a gadget block such as the mass or position arrays.
      //! Takes a lambda (or other function) which, given a mapper iterator, returns the data to be written for that
      //! particle. The writing proceeds in parallel using a memmap.
      template<typename WriteType>
      void saveGadgetBlock(particle::species forSpecies,
                           std::function<WriteType(const particle::mapper::MapperIterator<GridDataType> &)> getData) {

        size_t current_n = 0;
        size_t nTotalForThisBlock = 0;

        std::vector<decltype(writers[0].template getMemMapFortran<WriteType>(std::declval<size_t>()))> currentWriteBlocks;

        std::vector<size_t> nPerFile(nFiles,0);

        // What gadget particle types we are going to write?
        std::vector<unsigned int> particleTypes;
        for(unsigned int i=0; i<6; i++) {
          if(forSpecies == gadgetTypeToSpecies[i] || forSpecies == particle::species::all)
            particleTypes.push_back(i);
        }

        // Work out how many particles there will be in each file in this block
        for(auto pType: particleTypes) {
          for(unsigned int i=0; i<nFiles; i++) {
            nPerFile[i]+=this->nPartPerTypePerFile[i][pType];
            nTotalForThisBlock+=this->nPartPerTypePerFile[i][pType];
          }
        }

        if(forSpecies == particle::species::all)
          assert (nTotalForThisBlock==nTotal);

        for(int i=0; i<nFiles; i++)
          currentWriteBlocks.push_back(writers[i].template getMemMapFortran<WriteType>(nPerFile[i]));

        for (auto particle_type: particleTypes) {
          auto begin = mapper.beginParticleType(*generators[gadgetTypeToSpecies[particle_type]], particle_type);
          auto end = mapper.endParticleType(*generators[gadgetTypeToSpecies[particle_type]], particle_type);
          size_t nMax = end.getIndex() - begin.getIndex();

          current_n += begin.parallelIterate(
            [&](size_t n_offset, const particle::mapper::MapperIterator<GridDataType> &localIterator) {
              int fileNum = 0;
              size_t addr = n_offset + current_n;

              // Now figure out which file to dump this into.
              // The following approach looks slow (doing it for every particle!),
              // but code profiling suggests it's not a significant overhead.
              // Easier to do this than to try and have thread-local variables
              while(addr >= nPartPerFile[fileNum]) {
                addr-=nPartPerFile[fileNum];
                fileNum++;
              }

              currentWriteBlocks[fileNum][addr] = getData(localIterator);
            }, nMax);

        }

        assert(current_n == nTotalForThisBlock);

      }

      //! \brief Scan through data to obtain information on particles and their masses.
      void preScanForMassesAndParticleNumbers() {
        variableMass = false;
        masses = vector<InternalFloatType>(6, 0.0);
        nPartPerType = vector<size_t>(6, 0);
        nTotal = 0;

        logging::entry() << "Particles by gadget type:" << endl;

        for (unsigned int ptype = 0; ptype < 6; ++ptype) {
          InternalFloatType min_mass, max_mass;
          size_t n;
          getParticleInfo(min_mass, max_mass, n, ptype);
          if (n > 0) {
            if (min_mass != max_mass) {
              variableMass = true;
            }

            logging::entry() << "   Particle type " << ptype << ": " << n << " particles" << endl;

            nPartPerType[ptype] = size_t(n);
            masses[ptype] = min_mass;
            nTotal += n;
          }
        }

        if (variableMass) {
          logging::entry() << "Using variable-mass gadget format" << endl;
          for (unsigned int ptype = 0; ptype < 6; ++ptype) {
            masses[ptype] = 0;
          }
        } else {
          logging::entry() << "Using fixed-mass gadget format" << endl;
        }

        for(int i=0; i<nFiles; i++) {
          if (i == nFiles - 1) // last file
            nPartPerFile.push_back(nTotal - (nTotal / nFiles) * (nFiles - 1));
          else
            nPartPerFile.push_back(nTotal / nFiles);
        }

#ifdef DEBUG_INFO
        if(nFiles>1) {
          logging::entry() << "Particles per file: ";
          for (int i = 0; i < nFiles; i++)
            std::cerr << nPartPerFile[i] << " ";
          std::cerr << endl;
        }
#endif

      }

      //! \brief Extract the minimum mass, maximum mass, and number of a particle species
      void getParticleInfo(InternalFloatType &min_mass, InternalFloatType &max_mass,
                           size_t &num, unsigned int particle_type) {

        min_mass = std::numeric_limits<InternalFloatType>::max();
        max_mass = 0;
        num = 0;

        InternalFloatType mass;
        mapper.iterateParticlesOfType(*generators[gadgetTypeToSpecies[particle_type]], particle_type,
                                      [&](const auto &i) {
                                        mass = i.getMass(); // sometimes can be MUCH faster than getParticle
                                        if (min_mass > mass) min_mass = mass;
                                        if (max_mass < mass) max_mass = mass;
                                        num++;
                                      });

      }

      //! \brief Output the gadget3 or gadget2 header:
      void writeHeader() {
        std::vector<size_t> nPartPerTypeThisFile(6, 0);
        std::vector<size_t> nPartRemainingPerType = nPartPerType;
        size_t offset = 0;

        for(int i=0; i<nFiles; i++) {

          size_t nPartRemainingInFile = nPartPerFile[i];

          // now figure out what gadget particle types are going to appear in this particular file
          for(int partType = 0; partType<6 ; partType++) {
            if(nPartRemainingPerType[partType] < nPartRemainingInFile)
              nPartPerTypeThisFile[partType] = nPartRemainingPerType[partType];
            else
              nPartPerTypeThisFile[partType] = nPartRemainingInFile;

            nPartRemainingInFile-=nPartPerTypeThisFile[partType];
            nPartRemainingPerType[partType]-=nPartPerTypeThisFile[partType];
          }

#ifdef DEBUG_INFO
          logging::entry(logging::debug) << "i=" << i << " nPartPerTypeThisFile = ";
          for(int partType=0; partType<6; partType++) {
            std::cerr << nPartPerTypeThisFile[partType] << " ";
          }
          std::cerr << std::endl;
#endif


          assert(nPartRemainingInFile == 0); // ensure we have assigned all the particles to a type

          if (gadgetVersion == 3) {
            writers[i].writeFortran(createGadget3Header<OutputFloatType>(masses, nPartPerTypeThisFile, nPartPerType, nFiles,
                                                                         boxLength, cosmology));
          } else if (gadgetVersion == 2) {
            writers[i].writeFortran(createGadget2Header<OutputFloatType>(masses, nPartPerTypeThisFile, nPartPerType, nFiles,
                                                                         boxLength, cosmology));
          } else {
            throw std::runtime_error("Unknown gadget format");
          }

          this->nPartPerTypePerFile.push_back(nPartPerTypeThisFile);
        }

        // Cross-check that the various information we now hold on numbers of particles is self-consistent
        for(int partType=0; partType<6; partType++) {
          size_t numPerType = 0;
          for(int fileNum=0; fileNum<nFiles; fileNum++)
            numPerType+=nPartPerTypePerFile[fileNum][partType];
          assert(numPerType == this->nPartPerType[partType]);
        }
        for(int fileNum=0; fileNum<nFiles; fileNum++) {
          size_t numPerFile = 0;
          for(int partType=0; partType<6; partType++)
            numPerFile+=nPartPerTypePerFile[fileNum][partType];
          assert(numPerFile == this->nPartPerFile[fileNum]);
        }

      }

    public:
      /*! \brief Constructor

          \param boxLength - size of simulation box in Mpc/h.
          \param mapper - particle mapper used to create output particles.
          \param generators_ - vector of particle generators for each species.
          \param cosmology - struct containing cosmological parameters.
          \param gadgetVersion - 2 for gadget2 format, 3 for gadget3 format.
      */
      GadgetOutput(double boxLength,
                   particle::mapper::ParticleMapper<GridDataType> &mapper,
                   const particle::SpeciesToGeneratorMap<GridDataType> &generators_,
                   const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
                   int gadgetVersion, int numFiles) :
        mapper(mapper), generators(generators_), cosmology(cosmology), boxLength(boxLength),
        gadgetVersion(gadgetVersion), nFiles(numFiles) {
      }

      //! \brief Operation to save gadget particles
      void operator()(const std::string &name) {

        preScanForMassesAndParticleNumbers();

        if(nFiles==1) {
          writers.push_back(tools::MemMapFileWriter(name + std::to_string(gadgetVersion)));
        } else {
          for(int i=0; i<nFiles; i++) {
            writers.push_back(tools::MemMapFileWriter(name + std::to_string(gadgetVersion) + "." + std::to_string(i)));
          }
        }

        writeHeader();

        // positions
        saveGadgetBlock<Coordinate<OutputFloatType>>(
           particle::species::all,
          [](auto &localIterator) {
            auto particle = localIterator.getParticle();
            return Coordinate<OutputFloatType>(particle.pos);
          });

        // velocities
        saveGadgetBlock<Coordinate<OutputFloatType>>(
          particle::species::all,
          [](auto &localIterator) {
            auto particle = localIterator.getParticle();
            return Coordinate<OutputFloatType>(particle.vel);
          });

        // IDs
        saveGadgetBlock<unsigned long>(
          particle::species::all,
          [](auto &localIterator) {
            return localIterator.getIndex();
          });


        if (variableMass) {
          saveGadgetBlock<OutputFloatType>(
            particle::species::all,
            [](auto &localIterator) {
              return localIterator.getMass();
            });
        }

        if(cosmology.OmegaBaryons0>0) {
          // dummy internal energy
          saveGadgetBlock<OutputFloatType>(
            particle::species::baryon,
            [](auto &localIterator) -> OutputFloatType {
              return 0;
            });
        }

      }


    };


    //! \brief Creates GadgetOutput class and calls its save function.
    /*!
    \param name - name of output file
    \param Boxlength - simulation size in Mpc/h
    \param mapper - particle mapper used to link particles to grid locations
    \param generators - particles generators for each particle species (vector)
    \param cosmology - cosmological parameters
    \param gadgetformat - 2 or 3, gives type of gadget output (gadget2 or gadget3)
    */
    template<typename OutputFloatType, typename GridDataType>
    void save(const std::string &name, double Boxlength,
              particle::mapper::ParticleMapper<GridDataType> &mapper,
              particle::SpeciesToGeneratorMap<GridDataType> &generators,
              const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
              int gadgetformat, int nFiles) {

      GadgetOutput<GridDataType, OutputFloatType> output(Boxlength, mapper, generators, cosmology, gadgetformat, nFiles);
      output(name);

    }

  }
}


#endif
