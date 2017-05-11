#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <list>
#include "src/simulation/multilevelcontext/multilevelcontext.hpp"
#include "src/cosmology/parameters.hpp"
#include <src/simulation/particles/multilevelgenerator.hpp>

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

    template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
    class GraficOutput {
    protected:
      std::string outputFilename;
      multilevelcontext::MultiLevelContextInformation<DataType> context;
      std::shared_ptr<particle::AbstractMultiLevelParticleGenerator<DataType>> generator;
      const cosmology::CosmologicalParameters<T> &cosmology;

      T lengthFactor;
      T velFactor;
      size_t iordOffset;

    public:
      GraficOutput(const std::string &fname,
                   multilevelcontext::MultiLevelContextInformation<DataType> &levelContext,
                   particle::AbstractMultiLevelParticleGenerator<DataType> &particleGenerator,
                   const cosmology::CosmologicalParameters<T> &cosmology) :
        outputFilename(fname),
        generator(particleGenerator.shared_from_this()),
        cosmology(cosmology) {
        levelContext.copyContextWithIntermediateResolutionGrids(context);
        lengthFactor = 1. / cosmology.hubble; // Gadget Mpc a h^-1 -> GrafIC file Mpc a
        velFactor = std::pow(cosmology.scalefactor, 0.5f); // Gadget km s^-1 a^1/2 -> GrafIC km s^-1
      }

      void write() {
        iordOffset = 0;
        context.forEachLevel([&](const grids::Grid<T> &targetGrid) {
          writeGrid(targetGrid);
        });
      }

    protected:

      void writeGrid(const grids::Grid<T> &targetGrid) {
        auto evaluator = generator->makeEvaluatorForGrid(targetGrid);

        const grids::Grid<T> &baseGrid = context.getGridForLevel(0);
        size_t effective_size = tools::getRatioAndAssertPositiveInteger(baseGrid.dx * baseGrid.size, targetGrid.dx);
        tools::progress::ProgressBar pb("write grid " + std::to_string(effective_size), targetGrid.size);

        std::string thisGridFilename = outputFilename + "_" + std::to_string(effective_size);
        mkdir(thisGridFilename.c_str(), 0777);

        auto filenames = {"ic_velcx", "ic_velcy", "ic_velcz", "ic_poscx", "ic_poscy", "ic_poscz", "ic_particle_ids"};

        std::vector<size_t> block_lengths = {sizeof(float) * targetGrid.size2, sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2, sizeof(float) * targetGrid.size2,
                                             sizeof(float) * targetGrid.size2,
                                             sizeof(size_t) * targetGrid.size2};

        std::vector<std::ofstream> files;


        for (auto filename_i: filenames) {
          files.emplace_back(thisGridFilename + "/" + filename_i, std::ios::binary);
          writeHeaderForGrid(files.back(), targetGrid);
        }




        for (size_t i_z = 0; i_z < targetGrid.size; ++i_z) {
          pb.tick();
          writeBlockHeaderFooter(block_lengths, files);
          for (size_t i_y = 0; i_y < targetGrid.size; ++i_y) {
            for (size_t i_x = 0; i_x < targetGrid.size; ++i_x) {
              size_t i = targetGrid.getCellIndexNoWrap(i_x, i_y, i_z);
              size_t global_index = i + iordOffset;
              auto particle = evaluator->getParticleNoOffset(i);

              Coordinate<float> velScaled(particle.vel * velFactor);
              Coordinate<float> posScaled(particle.pos * lengthFactor);

              // Eek. The following code is horrible. Is there a way to make it neater?
              files[0].write((char *) (&velScaled.x), sizeof(float));
              files[1].write((char *) (&velScaled.y), sizeof(float));
              files[2].write((char *) (&velScaled.z), sizeof(float));
              files[3].write((char *) (&posScaled.x), sizeof(float));
              files[4].write((char *) (&posScaled.y), sizeof(float));
              files[5].write((char *) (&posScaled.z), sizeof(float));
              files[6].write((char *) (&global_index), sizeof(size_t));

            }
          }
          writeBlockHeaderFooter(block_lengths, files);
        }


        iordOffset += targetGrid.size3;

      }

      void writeBlockHeaderFooter(const vector<size_t> &block_lengths, vector<ofstream> &files) const {
        assert(block_lengths.size() == files.size());
        for (size_t i = 0; i < block_lengths.size(); ++i) {
          int block_length_as_integer = int(block_lengths[i]);
          files[i].write((char *) &block_length_as_integer, sizeof(int));
        }
      }

      void writeHeaderForGrid(std::ofstream &outFile, const grids::Grid<T> &targetGrid) {
        io_header_grafic header = getHeaderForGrid(targetGrid);
        int header_length = static_cast<int>(sizeof(io_header_grafic));
        outFile.write((char *) (&header_length), sizeof(int));
        outFile.write((char *) (&header), sizeof(io_header_grafic));
        outFile.write((char *) (&header_length), sizeof(int));
      }

      io_header_grafic getHeaderForGrid(const grids::Grid<T> &targetGrid) const {
        io_header_grafic header;
        header.nx = header.ny = header.nz = targetGrid.size;
        header.dx = targetGrid.dx * lengthFactor;
        header.xOffset = targetGrid.offsetLower.x * lengthFactor;
        header.yOffset = targetGrid.offsetLower.y * lengthFactor;
        header.zOffset = targetGrid.offsetLower.z * lengthFactor;
        header.scalefactor = cosmology.scalefactor;
        header.omegaM = cosmology.OmegaM0;
        header.omegaL = cosmology.OmegaLambda0;
        header.h0 = cosmology.hubble * 100;
        return header;
      }

    };

    template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
    void save(const std::string &filename,
              particle::AbstractMultiLevelParticleGenerator<DataType> &generator,
              multilevelcontext::MultiLevelContextInformation<DataType> &context,
              const cosmology::CosmologicalParameters<T> &cosmology) {
      GraficOutput<DataType> output(filename, context,
                                    generator, cosmology);
      output.write();
    }

  }
}
