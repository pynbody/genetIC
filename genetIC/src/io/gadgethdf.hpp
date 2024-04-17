#ifndef IC_GADGETHDF_HPP
#define IC_GADGETHDF_HPP

#include "gadget.hpp"
#include "H5Cpp.h"

namespace io {
  /*!
  \namespace io::gadget
  \brief Classes related to outputting particle data in the gadget 2 and gadget 3 formats.
*/
  namespace gadgethdf {


    using std::cerr;
    using std::endl;
    using std::vector;

    template<typename WriteType>
    constexpr const H5::PredType & getH5Type();

    template<>
    constexpr const H5::PredType & getH5Type<float>() {
      return H5::PredType::NATIVE_FLOAT;
    }

    template<>
    constexpr const H5::PredType & getH5Type<double>() {
      return H5::PredType::NATIVE_DOUBLE;
    }

    template<>
    constexpr const H5::PredType & getH5Type<unsigned long>() {
      return H5::PredType::NATIVE_ULONG;
    }

    template<>
    constexpr const H5::PredType & getH5Type<bool>() {
      return H5::PredType::NATIVE_HBOOL;
    }

    template<>
    constexpr const H5::PredType & getH5Type<unsigned int>() {
      return H5::PredType::NATIVE_UINT;
    }

    template<>
    constexpr const H5::PredType & getH5Type<Coordinate<double>>() {
      return getH5Type<double>();
    }

    template<>
    constexpr const H5::PredType & getH5Type<Coordinate<float>>() {
      return getH5Type<float>();
    }


    /*! \class GadgetOutput
    \brief Class to handle output to gadget files.
    */
    template<typename GridDataType, typename OutputFloatType>
    class GadgetHDFOutput : public gadget::GadgetOutput<GridDataType, OutputFloatType> {
    protected:
      using InternalFloatType = tools::datatypes::strip_complex<GridDataType>;
      std::vector<H5::H5File> h5Files;

      //! \brief Save a gadget block such as the mass or position arrays.
      //! Takes a lambda (or other function) which, given a mapper iterator, returns the data to be written for that
      //! particle. The writing proceeds in parallel using a memmap.
      template<typename WriteType>
      void saveGadgetHDFBlock(particle::species forSpecies,
                           std::function<WriteType(const particle::mapper::MapperIterator<GridDataType> &)> getData) {

      }

    public:

      GadgetHDFOutput(double boxLength,
                   particle::mapper::ParticleMapper<GridDataType> &mapper,
                   const particle::SpeciesToGeneratorMap<GridDataType> &generators_,
                   const cosmology::CosmologicalParameters<tools::datatypes::strip_complex<GridDataType>> &cosmology,
                   int numFiles) :
                   gadget::GadgetOutput<GridDataType, OutputFloatType>(boxLength, mapper, generators_,
                                                                       cosmology, 0, numFiles) { }



      template<typename WriteType, typename UnderlyingType=typename strip_coordinate<WriteType>::type>
      H5::DataSet createDataSet(H5::H5Location & location, std::string name, size_t nTotal) {
        const bool three_dimensional = ~std::is_same<UnderlyingType, WriteType>::value;

        if constexpr(three_dimensional) {
          static_assert(sizeof(Coordinate<UnderlyingType>) == 3*sizeof(UnderlyingType));
        }

        int rank;
        hsize_t dims_3d[2] = {nTotal, 3};
        hsize_t dims_1d[1] = {nTotal};
        hsize_t *dims;

        if(three_dimensional) {
          rank = 2;
          dims = dims_3d;
        } else {
          rank = 1;
          dims = dims_1d;
        }

        H5::DataSpace dataspace(rank, dims);
        H5::DataSet dataset = location.createDataSet(name, getH5Type<UnderlyingType>(), dataspace);
        return dataset;
      }


      template<typename WriteType>
      void saveBlock(particle::species forSpecies,
                     std::function<WriteType(const particle::mapper::MapperIterator<GridDataType> &)> getData,
                     std::string name) {

        size_t current_n = 0;
        size_t nTotalForThisBlock = 0;

        for(int particleType = 0; particleType < 6; particleType++) {
          if (forSpecies != this->gadgetTypeToSpecies[particleType] && forSpecies != particle::species::all)
            continue;

          std::vector<H5::Group> groups;
          std::vector<H5::DataSet> datasets;
          std::vector<WriteType> buffer;

          size_t nPartThisType = 0;

          // create or open particle type groups
          for (int i=0; i<this->nFiles; i++) {
            H5::H5File &file = h5Files[i];
            std::string groupName = "/PartType"+std::to_string(particleType);
            size_t nPart = this->nPartPerTypePerFile[i][particleType];
            try {
              groups.push_back(file.openGroup(groupName));
            } catch (H5::Exception &e) {
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
            datasets[i].write(&buffer[offset], getH5Type<WriteType>());
            offset+=this->nPartPerTypePerFile[i][particleType];
          }
        }

      }

      template<typename T>
      void writeHdfAttributeArray(H5::Group & group, std::string name, const std::vector<T> & value) {
        hsize_t size[1] = {value.size()};
        auto & type = getH5Type<T>();
        H5::Attribute attribute = group.createAttribute(name, type, H5::DataSpace(1, size));
        attribute.write(type, &value[0]);
      }

      template<typename T>
      void writeHdfAttribute(H5::Group & group, std::string name, const T & value) {
        auto & type = getH5Type<T>();
        H5::Attribute attribute = group.createAttribute(name, type, H5::DataSpace(H5S_SCALAR));
        attribute.write(type, &value);
      }

      virtual void writeHeaderOneFile(size_t fileNumber, std::vector<size_t> nPartPerTypeThisFile) {
        bool baryonic = (this->cosmology.OmegaBaryons0>0);
        H5::Group headerGroup = h5Files[fileNumber].createGroup("/Header");
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

        H5::Exception::dontPrint();

        this->preScanForMassesAndParticleNumbers();

        if(this->nFiles==1) {
          h5Files.push_back(H5::H5File(name + ".hdf5", H5F_ACC_TRUNC));
        } else {
          for(int i=0; i<this->nFiles; i++) {
            h5Files.push_back(H5::H5File(name + "." + std::to_string(i)+ ".hdf5", H5F_ACC_TRUNC));
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

      GadgetHDFOutput<GridDataType, OutputFloatType> output(Boxlength, mapper, generators, cosmology, nFiles);
      output(name);

    }

  }
}


#endif
