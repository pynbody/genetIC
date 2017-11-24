#ifndef IC_MASK_HPP
#define IC_MASK_HPP

#include "src/simulation/multilevelcontext/multilevelcontext.hpp"

namespace multilevelcontext {

  //! Abstract class to generate masks through the multilevel hierarchy
  template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
  class AbstractBaseMask {
  public:

    explicit AbstractBaseMask(MultiLevelContextInformation<DataType>* multilevelcontext_):
        multilevelcontext(multilevelcontext_), flaggedIdsAtEachLevel(multilevelcontext->getNumLevels()){}

    virtual T const isInMask(size_t level, size_t cellindex) = 0;

    virtual void ensureFlaggedVolumeIsContinuous() = 0;

    virtual void calculateMask() = 0;

  protected:
    MultiLevelContextInformation<DataType>* multilevelcontext;
    std::vector<std::vector<size_t>> flaggedIdsAtEachLevel;

    virtual void generateFlagsHierarchy() = 0;

  };

  //! Generate masks based on cells that have been flagged at each level.
  //! Useful for generating ic_refmap/ic_pvar files for RAMSES for example
  template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
  class Mask: public AbstractBaseMask<DataType, T> {

  public:
    explicit Mask(MultiLevelContextInformation<DataType>* multilevelcontext_):
        AbstractBaseMask<DataType, T> (multilevelcontext_){};

    /*!
     * @return 1.0 if cell with index is in the mask in this level, 0.0 else.
     */
    T const isInMask(size_t level, size_t cellindex) override{

      if(this->flaggedIdsAtEachLevel[level].size() == 0){
        // No flagged cells at all on this level, refine everywhere
        return T(1.0);
      } else if((std::binary_search(this->flaggedIdsAtEachLevel[level].begin(),
                                    this->flaggedIdsAtEachLevel[level].end(), cellindex))){
        // Cell is among flagged cells, refine it
        return T(1.0);
      }
      return T(0.0);

    }

    void calculateMask() override{
      generateFlagsHierarchy();
      ensureFlaggedVolumeIsContinuous();
      std::cerr << this->flaggedIdsAtEachLevel[0].size() << std::endl;
      std::cerr << this->flaggedIdsAtEachLevel[1].size() << std::endl;
      std::cerr << this->flaggedIdsAtEachLevel[2].size() << std::endl;
      for (auto i : this->flaggedIdsAtEachLevel[0]) {
        std::cerr << i << std::endl;
        std::cerr << this->multilevelcontext->getGridForLevel(0).getCellCoordinate(i) << std::endl;
      }
    }

    void ensureFlaggedVolumeIsContinuous() override{
      //TODO Add Lagrangian volume definition if it bothers you
    }

    void recalculateWithNewContext(MultiLevelContextInformation<DataType>* multilevelcontext_){
      this->multilevelcontext = multilevelcontext_;
      this->flaggedIdsAtEachLevel.clear();
      this->flaggedIdsAtEachLevel = std::vector<std::vector<size_t>>(this->multilevelcontext->getNumLevels());

      calculateMask();
    }

  protected:
    //! Calculates the flagged cells on all levels
    /*!
     * For each level, imports the flagged cells of the level if there are any. If not, e.g.
     * for virtual grids, downgrades the higher resolution mask on the virtual grid.
     */
    void generateFlagsHierarchy() override {
      size_t deepestFlaggedLevel;

      try {
        deepestFlaggedLevel = this->multilevelcontext->deepestLevelwithFlaggedCells();
      } catch (std::runtime_error& e){
        std::cerr << "WARNING No flagged particles were found on any level. Mask generation aborted" <<std::endl;
        return;
      }

      generateHierarchyAboveLevelInclusive(deepestFlaggedLevel);
      generateHierarchyBelowLevelExclusive(deepestFlaggedLevel);
    }

    void generateHierarchyAboveLevelInclusive(size_t deepestLevel){
      for (int level = deepestLevel; level >= 0; level--) {

        if (this->multilevelcontext->getGridForLevel(level).hasFlaggedCells()) {
          this->multilevelcontext->getGridForLevel(level).getFlaggedCells(this->flaggedIdsAtEachLevel[level]);

        } else {
          // Generate flags on intermediate levels that might not have some
          for (size_t i : this->flaggedIdsAtEachLevel[level + 1]) {
            this->flaggedIdsAtEachLevel[level].push_back(
                this->multilevelcontext->getIndexOfCellOnOtherLevel(level + 1, level, i));
          }
          // The above procedure might not be ordered and will create duplicates, get rid of them.
          sortAndEraseDuplicate(level);
        }
      }
    }

    void generateHierarchyBelowLevelExclusive(size_t coarsestLevel){

      // Do all level between this level and the finest
      for(size_t level=coarsestLevel + 1; level < this->multilevelcontext->getNumLevels() - 1; level++){

        for(size_t i = 0; i < this->multilevelcontext->getGridForLevel(level).size3; i++){
          size_t aboveindex = this->multilevelcontext->getIndexOfCellOnOtherLevel(level, level - 1, i);
          if(this->isMasked(aboveindex, level - 1)){
            this->flaggedIdsAtEachLevel[level].push_back(i);
          }
        }
      }
    }

    void sortAndEraseDuplicate(size_t level){
      tools::sortAndEraseDuplicate(this->flaggedIdsAtEachLevel[level]);
//      std::sort(this->flaggedIdsAtEachLevel[level].begin(), this->flaggedIdsAtEachLevel[level].end());
//      this->flaggedIdsAtEachLevel[level].erase(std::unique(
//          this->flaggedIdsAtEachLevel[level].begin(),
//          this->flaggedIdsAtEachLevel[level].end()), this->flaggedIdsAtEachLevel[level].end());
    }

    bool isMasked(size_t id, size_t level){
      return std::binary_search(this->flaggedIdsAtEachLevel[level].begin(),
                                this->flaggedIdsAtEachLevel[level].end(), id);
    }

  public:
    //! Generate a full multilevel field storing the mask information.
    /*! Increases memory requirement for the code but useful for debugging.
     */
    fields::MultiLevelField<DataType>* convertToField(){
      std::vector<std::shared_ptr<fields::Field<DataType, T>>> fields;

      // Field full of zeros
      for (size_t level = 0; level < this->multilevelcontext->getNumLevels(); ++level) {
        fields.push_back(std::shared_ptr<fields::Field<DataType, T>>(
            new fields::Field<DataType, T>(this->multilevelcontext->getGridForLevel(level))));
      }

      auto maskfield = new fields::MultiLevelField<DataType>(*(this->multilevelcontext), fields);

      // Field with mask information
      for (size_t level = 0; level < this->multilevelcontext->getNumLevels(); ++level) {
        for (size_t i=0; i< this->multilevelcontext->getGridForLevel(level).size3; i++){
          maskfield->getFieldForLevel(level).getDataVector()[i] = isInMask(level, i);
        }
      }

      return maskfield;
    }
  };
}



#endif
