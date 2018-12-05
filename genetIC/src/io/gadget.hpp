#ifndef IC_GADGET_HPP
#define IC_GADGET_HPP

#include "src/tools/memmap.hpp"
#include "src/tools/data_types/float_types.hpp"
#include "src/io.hpp"
#include <vector>


namespace io {
  namespace gadget {
    using std::cerr;
    using std::endl;
    using std::vector;


    struct io_header_2 //header for gadget2
    {
      int npart[6];
      double mass[6];
      double time;
      double redshift;
      int flag_sfr;
      int flag_feedback;
      int nPartTotal[6];
      /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the nPartTotal[2] stores the result of a division of the particle number by 2^32, while nPartTotal[1] holds the remainder. */
      int flag_cooling;
      int num_files;
      double BoxSize;
      double Omega0;
      double OmegaLambda;
      double HubbleParam;
      char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];  /* fills to 256 Bytes */
    };


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

/*!< holds header for snapshot files */



    template<typename OutputFloatType, typename InternalFloatType>
    io_header_2 createGadget2Header(vector<InternalFloatType> masses, vector<long> npart, double Boxlength,
                                    const cosmology::CosmologicalParameters<InternalFloatType> &cosmology) {
      io_header_2 header2;
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
      header2.nPartTotal[0] = npart[0];
      header2.nPartTotal[1] = npart[1];
      header2.nPartTotal[2] = npart[2];
      header2.nPartTotal[3] = npart[3];
      header2.nPartTotal[4] = npart[4];
      header2.nPartTotal[5] = npart[5];
      header2.flag_cooling = 0;
      header2.num_files = 1;
      header2.BoxSize = Boxlength;
      header2.Omega0 = cosmology.OmegaM0;
      header2.OmegaLambda = cosmology.OmegaLambda0;
      header2.HubbleParam = cosmology.hubble;

      if (npart[0] > 0) { //options for baryons
        header2.flag_sfr = 1;
        header2.flag_feedback = 1;
        header2.flag_cooling = 1;
      }


      return header2;
    }

    template<typename OutputFloatType, typename InternalFloatType>
    io_header_3 createGadget3Header(vector<InternalFloatType> masses, vector<long> npart, double Boxlength,
                                    const cosmology::CosmologicalParameters<InternalFloatType> &cosmology) {
      io_header_3 header3;
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
      header3.nPartTotal[0] = (unsigned int) (npart[0]);
      header3.nPartTotal[1] = (unsigned int) (npart[1]);
      header3.nPartTotal[2] = (unsigned int) (npart[2]);
      header3.nPartTotal[3] = (unsigned int) (npart[3]);
      header3.nPartTotal[4] = (unsigned int) (npart[4]);
      header3.nPartTotal[5] = (unsigned int) (npart[5]);
      header3.flag_cooling = 0;
      header3.num_files = 1;
      header3.BoxSize = Boxlength;
      header3.Omega0 = cosmology.OmegaM0;
      header3.OmegaLambda = cosmology.OmegaLambda0;
      header3.HubbleParam = cosmology.hubble;
      header3.flag_stellarage = 0;  /*!< flags whether the file contains formation times of star particles */
      header3.flag_metals = 0;    /*!< flags whether the file contains metallicity values for gas and star  particles */
      header3.nPartTotalHighWord[0] = (unsigned int) (npart[0] >> 32);
      header3.nPartTotalHighWord[1] = (unsigned int) (npart[1] >> 32); //copied from Gadget3
      header3.nPartTotalHighWord[2] = (unsigned int) (npart[2] >> 32);
      header3.nPartTotalHighWord[3] = (unsigned int) (npart[3] >> 32);
      header3.nPartTotalHighWord[4] = (unsigned int) (npart[4] >> 32);
      header3.nPartTotalHighWord[5] = (unsigned int) (npart[5] >> 32);
      header3.flag_entropy_instead_u = 0; /*!< flags that IC-file contains entropy instead of u */
      header3.flag_doubleprecision = tools::datatypes::floatinfo<OutputFloatType>::doubleprecision;
      header3.flag_ic_info = 1;
      header3.lpt_scalingfactor = 0.f; /*!dummy value since we never use ic_info!=1 */

      if (npart[0] > 0) { //options for baryons & special behavior
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

      particle::mapper::ParticleMapper<GridDataType> &mapper;
      std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>>> &generator;
      const cosmology::CosmologicalParameters<InternalFloatType> &cosmology;
      tools::MemMapFileWriter writer; // Used to write in the gadget format.
      size_t nTotal; // Total number of particles to output.
      double boxLength; // Size of simulation box.
      int gadgetVersion; // Which version of the gadget file to output. Allowed values 2 or 3.
      vector<InternalFloatType> masses; // Masses of particles if constant. Zero if variable.
      vector<long> npart; // Number of particles of each gadget type
      bool variableMass; // Stores whether we are using variable mass gadget particles


      // Map for the generators. Assume that the non-gas elements (gadget type 2+)
      // map to DM. Note that baryon generator is generator[1], DM is generator[0],
      // which is the opposite order to gadget.
      std::vector<size_t> gadgetTypeToFieldType {1,0,0,0,0,0};

      //! \brief Save a block of gadget particles
      template<typename WriteType>
      void saveGadgetBlock(std::function<WriteType(const particle::mapper::MapperIterator<GridDataType> &)> getData) {

        size_t current_n=0;

        auto currentWriteBlockC = writer.getMemMapFortran<WriteType>(nTotal);

        for(unsigned int particle_type=0; particle_type<6; particle_type++) {
          auto begin = mapper.beginParticleType(*generator[gadgetTypeToFieldType[particle_type]], particle_type);
          auto end = mapper.endParticleType(*generator[gadgetTypeToFieldType[particle_type]], particle_type);
          size_t nMax = end.getIndex()-begin.getIndex();

          current_n+=begin.parallelIterate([&](size_t n_offset, const particle::mapper::MapperIterator<GridDataType> &localIterator) {
            size_t addr = n_offset+current_n;
            currentWriteBlockC[addr] = getData(localIterator);
          }, nMax);

        }

        assert(current_n==nTotal);

      }

      //! \brief Scan through data to obtain information on particles and their masses.
      void preScanForMassesAndParticleNumbers() {
        variableMass = false;
        masses = vector<InternalFloatType>(6, 0.0);
        npart = vector<long>(6, 0);
        nTotal = 0;

        cerr << "Gadget output preliminary scan..." << endl;

        for (unsigned int ptype = 0; ptype < 6; ++ptype) {
          InternalFloatType min_mass, max_mass;
          size_t n;
          getParticleInfo(min_mass, max_mass, n, ptype);
          if (n > 0) {
            if (min_mass != max_mass) {
              variableMass = true;
            }

            cerr << "Particle type " << ptype << ": " << n << " particles" << endl;
            cerr << "Min and max particle mass : " << min_mass << " " << max_mass << endl;
            npart[ptype] = n;
            masses[ptype] = min_mass;
            nTotal += n;
          }
        }

        cerr << "" << endl;

        if (variableMass) {
          cerr << "Using variable-mass format" << endl;
          for (unsigned int ptype = 0; ptype < 6; ++ptype) {
            masses[ptype] = 0;
          }
        }
      }

        //! \brief Extract the minimum mass, maximum mass, and number of a particle species
      void getParticleInfo(InternalFloatType &min_mass, InternalFloatType &max_mass,
                           size_t &num, unsigned int particle_type) {

        min_mass = std::numeric_limits<InternalFloatType>::max();
        max_mass = 0;
        num = 0;

        InternalFloatType mass;
        mapper.iterateParticlesOfType(*generator[gadgetTypeToFieldType[particle_type]], particle_type, [&](const auto & i) {//CONFLICT_RESOLUTION
          mass = i.getMass(); // sometimes can be MUCH faster than getParticle
          if (min_mass > mass) min_mass = mass;
          if (max_mass < mass) max_mass = mass;
          num++;
        });

      }

      //! \brief Output the gadget3 or gadget2 header:
      void writeHeader() {
        if (gadgetVersion == 3) {
          writer.writeFortran(createGadget3Header<OutputFloatType>(masses, npart, boxLength, cosmology));
        } else if (gadgetVersion == 2) {
          writer.writeFortran(createGadget2Header<OutputFloatType>(masses, npart, boxLength, cosmology));
        } else {
          throw std::runtime_error("Unknown gadget format");
        }
      }

    public:
    //! \brief Constructor
      GadgetOutput(double boxLength,
                   particle::mapper::ParticleMapper<GridDataType> &mapper,
                   std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>>> &generator,
                   const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
                   int gadgetVersion) :
                   mapper(mapper), generator(generator), cosmology(cosmology), boxLength(boxLength),
                   gadgetVersion(gadgetVersion) {

                   //If not using gas, we should re-direct everything to DM:
                  if(generator.size() < 2)
                  {
                    for(auto& i : gadgetTypeToFieldType) {
                        i = 0;
                    }
                  }
       }

       //! \brief Operation to save gadget particles
       void operator()(const std::string &name) {

         preScanForMassesAndParticleNumbers();

         writer = tools::MemMapFileWriter(name + std::to_string(gadgetVersion));

         writeHeader();

         // positions
         saveGadgetBlock<Coordinate<OutputFloatType>>(
           [](auto &localIterator) {
             auto particle = localIterator.getParticle();
             return Coordinate<OutputFloatType>(particle.pos);
           });

         // velocities
         saveGadgetBlock<Coordinate<OutputFloatType>>(
           [](auto &localIterator) {
             auto particle = localIterator.getParticle();
             return Coordinate<OutputFloatType>(particle.vel);
           });

         // IDs
         saveGadgetBlock<long>(
           [](auto &localIterator) {
             return localIterator.getIndex();
           });


         if (variableMass) {
           saveGadgetBlock<OutputFloatType>(
             [](auto &localIterator) {
               return localIterator.getMass();
             });
         }
       }


    };


    //! \brief Creates GadgetOutput class and calls its save function.
    /*!
    \param name - name of output file
    \param Boxlength - simulation size in Mpc/h
    \param mapper - particle mapper used to link particles to grid locations
    \param generator - particles generators for each particle species (vector)
    \param cosmology - cosmological parameters
    \param gadgetformat - 2 or 3, gives type of gadget output (gadget2 or gadget3)
    */
    template<typename OutputFloatType, typename GridDataType>
    void save(const std::string &name, double Boxlength,
              particle::mapper::ParticleMapper<GridDataType> &mapper,
              //particle::AbstractMultiLevelParticleGenerator<GridDataType> &generator,//CONFLICT_RESOLUTION
              std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<GridDataType>>> &generator,//CONFLICT_RESOLUTION
              const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology, int gadgetformat) {

      GadgetOutput<GridDataType, OutputFloatType> output(Boxlength, mapper, generator, cosmology, gadgetformat);
      output(name);

    }

  }
}


#endif
