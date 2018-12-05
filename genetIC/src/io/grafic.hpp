#include <sys/stat.h>
#include <src/simulation/particles/multilevelgenerator.hpp>
#include <src/simulation/multilevelcontext/mask.hpp>
#include "src/tools/memmap.hpp"
#include <memory>
#include <vector>
#include <algorithm>
#include <string>

namespace io {
  namespace grafic {

    using std::ofstream;

    struct io_header_grafic {
      int nx, ny, nz;
      float dx;
      float xOffset, yOffset, zOffset;
      float scalefactor;
      float omegaM, omegaL, h0;
    } header_grafic;

    //! Export initial conditions in grafIC format, most likely for use with RAMSES.
    /**
     * WARNING: Grafic as described in Bertschinger 2001 uses "Mpc a" for header lengths and displacements, and
       "proper km s**-1" for velocities.
       However, RAMSES expects "Mpc a" for header, "Mpc a h^-1" for displacements and "proper km s**-1" for velocities,
       hence the need for the following three conversion factors.
       Beware of these units if using the displacements for other purposes than Ramses.
     */
    template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
    class GraficOutput {
    protected:
      std::string outputFilename;
      multilevelcontext::MultiLevelContextInformation<DataType> context;
      std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<DataType>>> generators;
      const cosmology::CosmologicalParameters<T> &cosmology;
      multilevelcontext::GraficMask<DataType,T>* mask;
      std::vector<std::shared_ptr<fields::OutputField<DataType>>> outputFields;


      T pvarValue;

      T lengthFactorHeader;
      T lengthFactorDisplacements;
      T velFactor;
      size_t iordOffset;

    public:
    //! \brief Constructor
      GraficOutput(const std::string &fname,
                   multilevelcontext::MultiLevelContextInformation<DataType> &levelContext,
                   std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<DataType>>> particleGenerators,
                   const cosmology::CosmologicalParameters<T> &cosmology,
                    const T pvarValue,
                   Coordinate<T> center,
                   size_t subsample,
                   size_t supersample,
                   std::vector<std::vector<size_t>>& input_mask,
                   std::vector<std::shared_ptr<fields::OutputField<DataType>>> outFields) :
          outputFilename(fname),
          cosmology(cosmology),
          pvarValue(pvarValue)
          {

        this->generators = particleGenerators;
        this->outputFields = outFields;

        levelContext.copyContextWithCenteredIntermediate(context, center, 2, subsample, supersample);

        mask = new multilevelcontext::GraficMask<DataType,T>(&context, input_mask);
        mask->calculateMask();

        lengthFactorHeader = 1. / cosmology.hubble; // Gadget Mpc a h^-1 -> GrafIC file Mpc a
        lengthFactorDisplacements = 1.;
        velFactor = std::pow(cosmology.scalefactor, 0.5f); // Gadget km s^-1 a^1/2 -> GrafIC km s^-1
      }

      //! \brief Output particles on all levels
      void write() {
        iordOffset = 0;

        for (size_t level = 0; level < this->context.getNumLevels(); ++level) {
            writeGrid(this->context.getGridForLevel(level), level);
        }
      }

    protected:

    //! \brief Output particles on a specified level
    /*!
    \param targetGrid - pointer to grid for specified level
    \param level - level to output.
    */
      void writeGrid(const grids::Grid<T> &targetGrid, size_t level) {

        auto evaluator_dm = generators[0]->makeParticleEvaluatorForGrid(targetGrid);
        auto overdensityFieldEvaluator = generators[0]->makeOverdensityEvaluatorForGrid(targetGrid);

	if(generators.size()>1) // use baryon field
	  overdensityFieldEvaluator = generators[1]->makeOverdensityEvaluatorForGrid(targetGrid);

        const grids::Grid<T> &baseGrid = context.getGridForLevel(0);
        size_t effective_size = tools::getRatioAndAssertPositiveInteger(baseGrid.cellSize * baseGrid.size,
                                                                        targetGrid.cellSize);
        tools::progress::ProgressBar pb("write grid " + std::to_string(effective_size), targetGrid.size);


        std::string thisGridFilename = outputFilename + "_" + std::to_string(effective_size);
        mkdir(thisGridFilename.c_str(), 0777);


        std::vector<std::string> filenames = {"ic_velcx", "ic_velcy", "ic_velcz", "ic_poscx", "ic_poscy",
                                              "ic_poscz",  "ic_deltab", "ic_refmap", "ic_pvar_00001", "ic_particle_ids"};

        std::vector<size_t> block_lengths = {sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(size_t) * targetGrid.size2};

        std::vector<tools::MemMapFileWriter> files;



        for (size_t i=0; i<filenames.size(); ++i) {
          auto filename_i = filenames[i];
          files.emplace_back(thisGridFilename + "/" + filename_i);
          writeHeaderForGrid(files.back(), targetGrid);
        }

        std::shared_ptr<fields::Field<DataType, T>> baryonFieldOnLevelPtr = nullptr;
        if(this->outputFields.size() > 1)
        {
            for(size_t i = 0; i < outputFields.size(); i++)
            {
                this->outputFields[i]->toReal(); // Field should be output in real space
            }

            baryonFieldOnLevelPtr = outputFields[1]->getFieldForLevel(level).shared_from_this();
        }

        for (size_t i_z = 0; i_z < targetGrid.size; ++i_z) {
          pb.tick();
          writeBlockHeaderFooter(block_lengths, files);

          std::vector<tools::MemMapRegion<float>> varMaps;
          for(int m=0; m<9; ++m)
            varMaps.push_back(files[m].getMemMap<float>(targetGrid.size2));

          tools::MemMapRegion<size_t> idMap = files[9].getMemMap<size_t>(targetGrid.size2);

#pragma omp parallel for
          for (size_t i_y = 0; i_y < targetGrid.size; ++i_y) {
            for (size_t i_x = 0; i_x < targetGrid.size; ++i_x) {
              size_t i = targetGrid.getIndexFromCoordinateNoWrap(i_x, i_y, i_z);
              size_t global_index = i + iordOffset;
              auto particle = evaluator_dm->getParticleNoOffset(i);

              Coordinate<float> velScaled(particle.vel * velFactor);
              Coordinate<float> posScaled(particle.pos * lengthFactorDisplacements);


              float deltab = (*overdensityFieldEvaluator)[i];

              // Detect whether we are using baryons:
              float mask = this->mask->isInMask(level, i);
              float pvar = pvarValue * mask;
              size_t file_index = i_y*targetGrid.size+i_x;


              varMaps[0][file_index] = velScaled.x;
              varMaps[1][file_index] = velScaled.y;
              varMaps[2][file_index] = velScaled.z;
              varMaps[3][file_index] = posScaled.x;
              varMaps[4][file_index] = posScaled.y;
              varMaps[5][file_index] = posScaled.z;
              varMaps[6][file_index] = deltab;
              varMaps[7][file_index] = mask;
              varMaps[8][file_index] = pvar;
              idMap[file_index] = global_index;

            }
          }
          writeBlockHeaderFooter(block_lengths, files);
        }
        iordOffset += targetGrid.size3;
      }

      //! \brief Output specified blocks of data to specified files.
      /*!
      \param block_lengths - lengths of blocks of data to output, for each file.
      \param files - files to output to.
      */
      void writeBlockHeaderFooter(const vector<size_t> &block_lengths, vector<tools::MemMapFileWriter> &files) const {
        assert(block_lengths.size() == files.size());
        for (size_t i = 0; i < block_lengths.size(); ++i) {
          files[i].write<int>(int(block_lengths[i]));
        }
      }

      //! \brief Output the header for a given level of the simulation.
      /*!
      \param file - file to output the header to
      \param targetGrid - pointer to the grid to output header for
      */
      void writeHeaderForGrid(tools::MemMapFileWriter & file, const grids::Grid<T> &targetGrid) {
        io_header_grafic header = getHeaderForGrid(targetGrid);
        int header_length = static_cast<int>(sizeof(io_header_grafic));
        file.write<int>(header_length);
        file.write<io_header_grafic>(header);
        file.write<int>(header_length);
      }

      //! \brief Returns the grafic header for the specified level of the simulation
      /*!
      \param targetGrid - pointer to the grid to retrieve the header for
      */
      io_header_grafic getHeaderForGrid(const grids::Grid<T> &targetGrid) const {
        io_header_grafic header;
        header.nx = header.ny = header.nz = targetGrid.size;
        header.dx = targetGrid.cellSize * lengthFactorHeader;
        header.xOffset = targetGrid.offsetLower.x * lengthFactorHeader;
        header.yOffset = targetGrid.offsetLower.y * lengthFactorHeader;
        header.zOffset = targetGrid.offsetLower.z * lengthFactorHeader;
        header.scalefactor = cosmology.scalefactor;
        header.omegaM = cosmology.OmegaM0;
        header.omegaL = cosmology.OmegaLambda0;
        header.h0 = cosmology.hubble * 100;
        return header;
      }

    };

    //! \brief Save particles in grafic format.
    /*!
    \param filename - name of output file
    \param generators - vector of particles generators for each type of particle (dark matter, baryons)
    \param context - multi-level context on which the overdensity fields are defined.
    \param cosmology - cosmological parameters
    \param pvarValue - passive variables for refinement masks
    \param center - centre point to use for centred grids
    \param subsample - subsample factor specified in paramter file
    \param supersample - supersample factor specified in paramter file
    \param input_mask - Grafic mask being used
    \param outputFields - Vector of references to underlying overdensity fields
    */
    template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
    void save(const std::string &filename,
              std::vector<std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<DataType>>> generators,
              multilevelcontext::MultiLevelContextInformation<DataType> &context,
              const cosmology::CosmologicalParameters<T> &cosmology,
              const T pvarValue, Coordinate<T> center, size_t subsample, size_t supersample,
              std::vector<std::vector<size_t>>& input_mask,
              std::vector<std::shared_ptr<fields::OutputField<DataType>>>& outputFields) {
      GraficOutput<DataType> output(filename, context, generators,
                                    cosmology, pvarValue, center, subsample, supersample, input_mask,outputFields);
      output.write();
    }

  }
}
