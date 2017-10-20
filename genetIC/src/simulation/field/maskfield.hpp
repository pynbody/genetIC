#ifndef IC_MASKFIELD_HPP
#define IC_MASKFIELD_HPP


#include "src/simulation/field/multilevelfield.hpp"
namespace fields {

  template<typename DataType>
  class AbstractMaskField : public MultiLevelField<DataType> {

  public:
    explicit AbstractMaskField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext)
        : MultiLevelField<DataType>(multiLevelContext) {}

    virtual void calculateMasksAllLevels() = 0;

  protected:
    using T = typename MultiLevelField<DataType>::T;

    //! From flagged cells, calculate refinement masks, i.e. 0 and 1 maps
    virtual void calculateMaskFinestLevel() = 0;

    //! Find center of volume defined by flagged cells
    virtual Coordinate<T> findVolumeCenter() = 0;

    //! Flag supplementary cells to have a continuous volume and/or define binding volume as sphere/cube
    virtual void flagBindingVolume() = 0;

  private:
    virtual Coordinate<T> calculateVolumeCenter(std::vector<size_t>) = 0;

  };

  template<typename DataType>
  class MaskField : public AbstractMaskField<DataType> {

  private:
    std::vector<size_t> flaggedcellsarray;

  public:
    explicit MaskField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext)
        : AbstractMaskField<DataType>(multiLevelContext){

      size_t finestlevel = this->multiLevelContext->getNumLevels() - 1;
      auto finestgrid = this->multiLevelContext->getGridForLevel(finestlevel);
      finestgrid.getFlaggedCells(flaggedcellsarray);
    }


    void calculateMasksAllLevels() override {};

  protected:
    using T = typename AbstractMaskField<DataType>::T;

    void calculateMaskFinestLevel() override {
      size_t finestlevel = this->multiLevelContext->getNumLevels() - 1;

      Coordinate<T> centre = this->findVolumeCenter();
      this->flagBindingVolume();

      // Zero the entire field
      for (size_t i = 0; i < this->multiLevelContext->getGridForLevel(finestlevel).size3; ++i) {
        this->fieldsOnLevels[finestlevel][i] = 0;
      }

      // Allow (=1) refinement on flagged cells
      for (size_t i = 0; i < flaggedcellsarray.size(); ++i) {
        this->fieldsOnLevels[finestlevel][flaggedcellsarray[i]] = 1;
      }
    }

    Coordinate<T> findVolumeCenter() override {

    };

    void flagBindingVolume() override {};

  private:
    Coordinate<T> calculateVolumeCenter(std::vector<size_t>) override{

    }
  };




//  template<typename DataType>
//  class ConvexHullMaskField : public MaskField<DataType> {
//
//  public:
//    explicit ConvexHullMaskField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext)
//        : MaskField<DataType>(multiLevelContext){}
//
//  };

  //! Not sure I actually need this since grid might implement most of this functionalities
  class AbstractVolumeCalculator {

  protected:
    virtual void findCenter() =0;

    virtual void extractBindingVolume() = 0;

  };

  class ConvexHullVolumeCalculator: AbstractVolumeCalculator{

  protected:
    void findCenter() override{}

    void extractBindingVolume() override{}
  };




}
#endif
