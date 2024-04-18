#ifndef IC_GADGETHDF_HPP
#define IC_GADGETHDF_HPP

#include "gadget.hpp"

// we use this third-party HDF C++ wrapper because it is more modern and doesn't require
// a binary C++ dylib. The C++ dylib causes major problems on macos trying to get everything
// to link correctly (e.g. GNU vs clang libstdc++ issues)
#include <highfive/highfive.hpp>

namespace io {
  /*!
  \namespace io::gadget
  \brief Classes related to outputting particle data in the gadget 2 and gadget 3 formats.
*/
  namespace gadgethdf {


    using std::cerr;
    using std::endl;
    using std::vector;


    /*! \class GadgetOutput
    \brief Class to handle output to gadget files.
    */
    template<typename GridDataType, typename OutputFloatType>
    class GadgetHDFOutput : public gadget::GadgetOutput<GridDataType, OutputFloatType> {
    protected:
      using InternalFloatType = tools::datatypes::strip_complex<GridDataType>;
      std::vector<HighFive::File> h5Files;

    public:

      GadgetHDFOutput(double boxLength,
                   particle::mapper::ParticleMapper<GridDataType> &mapper,
                   const particle::SpeciesToGeneratorMap<GridDataType> &generators_,
                   const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
                   int numFiles) :
                   gadget::GadgetOutput<GridDataType, OutputFloatType>(boxLength, mapper, generators_,
                                                                       cosmology, 0, numFiles) { }



      template<typename WriteType, typename UnderlyingType=typename strip_coordinate<WriteType>::type>
      HighFive::DataSet createDataSet(HighFive::Group & location, std::string name, size_t nTotal) {
        const bool three_dimensional = !std::is_same<UnderlyingType, WriteType>::value;

        std::unique_ptr<HighFive::DataSpace> ds;

        if constexpr(three_dimensional) {
          static_assert(sizeof(Coordinate<UnderlyingType>) == 3*sizeof(UnderlyingType));
          ds = std::make_unique<HighFive::DataSpace>(nTotal, 3);
        } else {
          ds = std::make_unique<HighFive::DataSpace>(nTotal);
        }

        HighFive::DataSet dataset = location.createDataSet<UnderlyingType>(name, *ds);
        return dataset;
      }


      template<typename WriteType>
      void saveBlock(particle::species forSpecies,
                     std::function<WriteType(const particle::mapper::MapperIterator<GridDataType> &)> getData,
                     std::string name) {

        using UnderlyingType = typename strip_coordinate<WriteType>::type;

        size_t current_n = 0;
        size_t nTotalForThisBlock = 0;

        for(int particleType = 0; particleType < 6; particleType++) {
          if (forSpecies != this->gadgetTypeToSpecies[particleType] && forSpecies != particle::species::all)
            continue;

          std::vector<HighFive::Group> groups;
          std::vector<HighFive::DataSet> datasets;
          std::vector<WriteType> buffer;

          size_t nPartThisType = 0;

          // create or open particle type groups
          for (int i=0; i<this->nFiles; i++) {
            HighFive::File &file = h5Files[i];
            std::string groupName = "/PartType"+std::to_string(particleType);
            size_t nPart = this->nPartPerTypePerFile[i][particleType];
            try {
              groups.push_back(file.getGroup(groupName));
            } catch (HighFive::Exception &e) {
              groups.push_back(file.createGroup(groupName));
            }
            datasets.push_back(createDataSet<WriteType>(groups.back(), name, nPart));
            nPartThisType += nPart;
          }

          assert(nPartThisType == this->nPartPerType[particleType]);

          buffer.resize(nPartThisType);

          auto begin = this->mapper.beginParticleType(*this->generators[this->gadgetTypeToSpecies[particleType]],
                                                      particleType);
          auto end = this->mapper.endParticleType(*this->generators[this->gadgetTypeToSpecies[particleType]],
                                                  particleType);
          begin.parallelIterate(
              [&](size_t n_offset, const particle::mapper::MapperIterator<GridDataType> &localIterator) {
                buffer[n_offset] = getData(localIterator);
              }, end - begin);

          size_t offset = 0;
          for (int i = 0; i < this->nFiles; i++) {
            datasets[i].write_raw<UnderlyingType>(reinterpret_cast<UnderlyingType*>(&buffer[offset]));
            offset+=this->nPartPerTypePerFile[i][particleType];
          }
        }

      }

      template<typename T>
      void writeHdfAttributeArray(HighFive::Group & group, std::string name, const std::vector<T> & value) {
        HighFive::Attribute attribute = group.createAttribute<T>(name, HighFive::DataSpace(value.size()));
        attribute.write( value);
      }

      template<typename T>
      void writeHdfAttribute(HighFive::Group & group, std::string name, const T & value) {
        HighFive::Attribute attribute = group.createAttribute<T>(name,
                   HighFive::DataSpace(HighFive::DataSpace::DataspaceType::dataspace_scalar));
        attribute.write<T>( value);
      }

      virtual void writeHeaderOneFile(size_t fileNumber, std::vector<size_t> nPartPerTypeThisFile) {
        bool baryonic = (this->cosmology.OmegaBaryons0>0);
        HighFive::Group headerGroup = h5Files[fileNumber].createGroup("/Header");
        writeHdfAttributeArray(headerGroup, "NumPart_ThisFile", nPartPerTypeThisFile);
        writeHdfAttributeArray(headerGroup, "NumPart_Total", this->nPartPerType);
        writeHdfAttributeArray(headerGroup, "MassTable", this->masses);
        writeHdfAttribute(headerGroup, "ExpansionFactor", this->cosmology.scalefactor);
        writeHdfAttribute(headerGroup, "Redshift", 1./this->cosmology.scalefactor-1.);
        writeHdfAttribute(headerGroup, "HubbleParam", this->cosmology.hubble);
        writeHdfAttribute(headerGroup, "Flag_Sfr", baryonic);
        writeHdfAttribute(headerGroup, "Flag_Cooling", baryonic);
        writeHdfAttribute(headerGroup, "Flag_Feedback", baryonic);
        writeHdfAttribute(headerGroup, "Flag_StellarAge", baryonic);
        writeHdfAttribute(headerGroup, "Flag_Metals", baryonic);
        writeHdfAttribute(headerGroup, "BoxSize", this->boxLength);
        writeHdfAttribute(headerGroup, "Omega0", this->cosmology.OmegaM0);
        writeHdfAttribute(headerGroup, "OmegaLambda", this->cosmology.OmegaLambda0);
        writeHdfAttribute(headerGroup, "Flag_DoublePrecision", std::is_same<OutputFloatType, double>::value);
        writeHdfAttribute(headerGroup, "NumFilesPerSnapshot", this->nFiles);



      }


      //! \brief Operation to save gadget particles
      virtual void operator()(const std::string &name) {

        this->preScanForMassesAndParticleNumbers();

        if(this->nFiles==1) {
          h5Files.push_back(HighFive::File(name + ".hdf5", HighFive::File::Truncate));
        } else {
          for(int i=0; i<this->nFiles; i++) {
            h5Files.push_back(HighFive::File(name + "." + std::to_string(i)+ ".hdf5", HighFive::File::Truncate));
          }
        }


        this->writeHeader();



        saveBlock<Coordinate<OutputFloatType>>(
           particle::species::all,
          [](auto &localIterator) {
            auto particle = localIterator.getParticle();
            return Coordinate<OutputFloatType>(particle.pos);
          }, "Coordinates");

        saveBlock<Coordinate<OutputFloatType>>(
          particle::species::all,
          [](auto &localIterator) {
            auto particle = localIterator.getParticle();
            return Coordinate<OutputFloatType>(particle.vel);
          }, "Velocities");

        saveBlock<unsigned long>(
          particle::species::all,
          [](auto &localIterator) {
            return localIterator.getIndex();
          }, "ParticleIDs");


        if (this->variableMass) {
          saveBlock<OutputFloatType>(
            particle::species::all,
            [](auto &localIterator) {
              return localIterator.getMass();
            }, "Masses");
        }




      }


    };


    //! \brief Generates gadget HDF output
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
      GadgetHDFOutput<GridDataType, OutputFloatType> output(Boxlength, mapper, generators, cosmology, nFiles);
      output(name);

    }

  }
}


#endif
