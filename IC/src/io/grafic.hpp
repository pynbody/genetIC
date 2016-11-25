//
// Created by Andrew Pontzen on 17/11/2016.
//

#include <fstream>
#include <cmath>
#include "../multilevelcontext.hpp"
#include "../cosmo.hpp"
#include <sys/stat.h>
#include <sys/types.h>

namespace io {
  namespace grafic {

    struct io_header_grafic {
      int nx, ny, nz;
      float dx;
      float xOffset, yOffset, zOffset;
      float scalefactor;
      float omegaM, omegaL, h0;
    } header_grafic;

    template<typename T>
    class GraficOutput {
    protected:
      std::string outputFilename;
      MultiLevelContextInformation<T> context;
      const CosmologicalParameters<T> &cosmology;

      T lengthFactor;
      T velFactor;
      size_t iordOffset;

    public:
      GraficOutput(const std::string & fname,
                   MultiLevelContextInformation<T> & levelContext,
                   const CosmologicalParameters<T> &cosmology):
        outputFilename(fname),
        cosmology(cosmology)
      {
        levelContext.copyContextWithIntermediateResolutionGrids(context);
        lengthFactor = 1./cosmology.hubble; // Gadget Mpc a h^-1 -> GrafIC file Mpc a
        velFactor = std::pow(cosmology.scalefactor, 0.5f); // Gadget km s^-1 a^1/2 -> GrafIC km s^-1
      }

      void write() {
        iordOffset=0;
        context.forEachLevel([&](const Grid<T> & targetGrid) {
          writeGrid(targetGrid);
        });
      }

    protected:
      static size_t getRatioAndAssertInteger(T a, T b, T tol=1e-8) {
        T ratio = a/b;
        size_t ratio_int = size_t(round(ratio));

        T residual = T(ratio_int)-ratio;
        assert(abs(residual)<tol);
        return ratio_int;
      }

      void writeGrid(const Grid<T> & targetGrid) {

        const Grid<T> & baseGrid = context.getGridForLevel(0);
        size_t effective_size =  getRatioAndAssertInteger(baseGrid.dx*baseGrid.size, targetGrid.dx);
	progress::ProgressBar pb("write grid "+std::to_string(effective_size), targetGrid.size3);

        std::string thisGridFilename = outputFilename+"_"+std::to_string(effective_size);
        mkdir(thisGridFilename.c_str(), 0777);

        std::ofstream outFileX(thisGridFilename+"/ic_velcx", ios::binary);
        std::ofstream outFileY(thisGridFilename+"/ic_velcy", ios::binary);
        std::ofstream outFileZ(thisGridFilename+"/ic_velcz", ios::binary);
        std::ofstream outFilePID(thisGridFilename+"/ic_particle_ids", ios::binary);

        writeHeaderForGrid(outFileX, targetGrid);
        writeHeaderForGrid(outFileY, targetGrid);
        writeHeaderForGrid(outFileZ, targetGrid);
        writeHeaderForGrid(outFilePID, targetGrid);

        int data_block_length = static_cast<int>(sizeof(float)*targetGrid.size2);
        int iord_block_length = static_cast<int>(sizeof(int)*targetGrid.size2);

        outFileX.write((char*)(&data_block_length), sizeof(int));
        outFileY.write((char*)(&data_block_length), sizeof(int));
        outFileZ.write((char*)(&data_block_length), sizeof(int));
        outFilePID.write((char*)(&iord_block_length), sizeof(int));

        for(size_t i=0; i<targetGrid.size3; ++i) {
	  pb.tick();
          size_t global_index = i + iordOffset;
          auto particle = targetGrid.getParticleNoOffset(i);
          Coordinate<float> velScaled = particle.vel*velFactor;

          outFileX.write((char*)(&velScaled.x), sizeof(float));
          outFileY.write((char*)(&velScaled.y), sizeof(float));
          outFileZ.write((char*)(&velScaled.z), sizeof(float));
          outFilePID.write((char*)(&global_index), sizeof(size_t));

          if(i%targetGrid.size2==targetGrid.size2-1 && i<targetGrid.size3-2) {
            // end of one "slice", start a new fortran block in the file
            outFileX.write((char*)(&data_block_length), sizeof(int));
            outFileY.write((char*)(&data_block_length), sizeof(int));
            outFileZ.write((char*)(&data_block_length), sizeof(int));
            outFilePID.write((char*)(&iord_block_length), sizeof(int));

            outFileX.write((char*)(&data_block_length), sizeof(int));
            outFileY.write((char*)(&data_block_length), sizeof(int));
            outFileZ.write((char*)(&data_block_length), sizeof(int));
            outFilePID.write((char*)(&iord_block_length), sizeof(int));
          }

        }

        outFileX.write((char*)(&data_block_length), sizeof(int));
        outFileY.write((char*)(&data_block_length), sizeof(int));
        outFileZ.write((char*)(&data_block_length), sizeof(int));
        outFilePID.write((char*)(&iord_block_length), sizeof(int));

        iordOffset+=targetGrid.size3;

      }

      void writeHeaderForGrid(std::ofstream &outFile, const Grid<T> &targetGrid) {
        io_header_grafic header = getHeaderForGrid(targetGrid);
        int header_length = static_cast<int>(sizeof(io_header_grafic));
        outFile.write((char*)(&header_length), sizeof(int));
        outFile.write((char*)(&header), sizeof(io_header_grafic));
        outFile.write((char*)(&header_length), sizeof(int));
      }

      io_header_grafic getHeaderForGrid(const Grid<T> &targetGrid) const {
        io_header_grafic header;
        header.nx = header.ny = header.nz = targetGrid.size;
        header.dx = targetGrid.dx * lengthFactor;
        header.xOffset = targetGrid.offsetLower.x * lengthFactor;
        header.yOffset = targetGrid.offsetLower.y * lengthFactor;
        header.zOffset = targetGrid.offsetLower.z * lengthFactor;
        header.scalefactor = cosmology.scalefactor;
        header.omegaM = cosmology.OmegaM0;
        header.omegaL = cosmology.OmegaLambda0;
        header.h0 = cosmology.hubble*100;
        return header;
      }

    };

    template<typename T>
    void save(const std::string & filename,  MultiLevelContextInformation<T> &context,
              const CosmologicalParameters<T> &cosmology) {
      GraficOutput<T> output(filename,context,cosmology);
      output.write();
    }

  }
}
