#ifndef IC_MASKFIELD_HPP
#define IC_MASKFIELD_HPP

#include "src/simulation/field/multilevelfield.hpp"


namespace fields {

  //! Defining multilevel masks, i.e. fields of 0 and 1 to track where something is allowed on a multilevel hierarchy
  template<typename DataType>
  class AbstractMaskField : public MultiLevelField<DataType> {

  public:

    //! Define a MultiLevelField full of zeroes
    explicit AbstractMaskField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext):
        MultiLevelField<DataType>(multiLevelContext){

      std::vector<std::shared_ptr<Field<DataType, T>>> fields;

      for (size_t level = 0; level < multiLevelContext.getNumLevels(); ++level) {
          fields.push_back(std::shared_ptr<Field<DataType, T>>(
                            new Field<DataType, T>(multiLevelContext.getGridForLevel(level))));
      }

      this->fieldsOnLevels = fields;
    }

    virtual void calculateMasksAllLevels() = 0;

  protected:
    using T = typename MultiLevelField<DataType>::T;

    virtual void calculateMaskFinestLevel() = 0;

    //! Find center of volume defined by flagged cells
    virtual Coordinate<T> findVolumeCenter() = 0;

    //! Flag supplementary cells to have a continuous volume and/or define binding volume as sphere/cube
    virtual void flagBindingVolume() = 0;

  private:
    virtual Coordinate<T> calculateVolumeCenter(std::vector<size_t>) = 0;

  };

  //! This implementation of mask is targeted for RAMSES refinement masks.
  template<typename DataType>
  class RAMSESMaskField : public AbstractMaskField<DataType> {

  private:
    std::vector<size_t> flaggedcellsarray;

  public:
    explicit RAMSESMaskField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext)
        : AbstractMaskField<DataType>(multiLevelContext){

      size_t finestlevel = this->multiLevelContext->getNumLevels() - 1;
      auto finestgrid = this->multiLevelContext->getGridForLevel(finestlevel);
      finestgrid.getFlaggedCells(flaggedcellsarray);
    }

    //!
    // TODO: Implementing for now the easy way
    void calculateMasksAllLevels() override {

      if(this->multiLevelContext->getNumLevels() == 1){
        calculateMaskFinestLevel();
      }
      else {
        calculateMaskFinestLevel();


        for (size_t level = 0; level < this->multiLevelContext->getNumLevels() - 1; ++level) {

          // Loop over all cells and see if it is contained in the next level grid
          // grid.containsPoint might be the useful method
          // This loop should be easily pragma openmped

          auto current_grid = this->multiLevelContext->getGridForLevel(level);
          auto next_level_grid = this->multiLevelContext->getGridForLevel(level + 1);

#pragma omp parallel for
          for (size_t i = 0; i < current_grid.size3; i++) {
            Coordinate<int> cell_coord = current_grid.getCellCoordinate(i);

            if (next_level_grid.containsCellWithCoordinate(cell_coord)) {
              this->getFieldForLevel(level).getDataVector()[i] = 1;
            } else {
              this->getFieldForLevel(level).getDataVector()[i] = 0;
            }
          }

        }
      }

    };

  protected:
    using T = typename AbstractMaskField<DataType>::T;

    //! Calculate finest level mask as flagged particles as this level
    void calculateMaskFinestLevel() override {
      size_t finestlevel = this->multiLevelContext->getNumLevels() - 1;

//      Coordinate<T> centre = this->findVolumeCenter();
      this->flagBindingVolume();

      // Allow (=1) refinement on flagged cells
      for (size_t i = 0; i < flaggedcellsarray.size(); ++i) {
        this->getFieldForLevel(finestlevel).getDataVector()[flaggedcellsarray[i]] = 1;
      }
    }

    //! Find the center defined by the array of flagged cells
    Coordinate<T> findVolumeCenter() override {
      return Coordinate<T>();

    };

    //! The volume defined by the array might have gaps in it. Deal with this.
    void flagBindingVolume() override {};

  private:
    Coordinate<T> calculateVolumeCenter(std::vector<size_t>) override{
      return Coordinate<T>();
    }
  };

  //! Not sure I actually need this since grid might implement most of this functionalities
  // ideally the definition of the volume would be an independent object
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
