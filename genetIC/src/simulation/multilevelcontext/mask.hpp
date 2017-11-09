#ifndef IC_MASK_HPP
#define IC_MASK_HPP

#include "src/simulation/multilevelcontext/multilevelcontext.hpp"
#include <cstdlib>
#include <vector>

namespace multilevelcontext {

//  template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
//  class MultiLevelContextInformation<DataType>;

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

  template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
  class RamsesMask: public AbstractBaseMask<DataType, T> {

  public:
    explicit RamsesMask(MultiLevelContextInformation<DataType>* multilevelcontext_):
        AbstractBaseMask<DataType, T> (multilevelcontext_){};

    //! Mask generation if this is useful for your application, e.g. Ramses
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
    void generateFlagsHierarchy() override {

      // Mask is generated assuming the deepest flagged cells are the interesting cells one is after
      size_t deepestFlaggedLevel = this->multilevelcontext->deepestLevelwithFlaggedCells();
      generateHierarchyAboveLevelInclusive(deepestFlaggedLevel);
      generateHierarchyBelowLevelExclusive(deepestFlaggedLevel);
    }

    void generateHierarchyAboveLevelInclusive(size_t deepestLevel){
      for (int level = deepestLevel; level >= 0; level--) {

        auto currentLevelGrid = this->multilevelcontext->getGridForLevel(level);
        if (currentLevelGrid.hasFlaggedCells()) {
          currentLevelGrid.getFlaggedCells(this->flaggedIdsAtEachLevel[level]);

        } else {
          // Generate flags on intermediate levels that might not have some
          for (size_t i : this->flaggedIdsAtEachLevel[level + 1]) {
            this->flaggedIdsAtEachLevel[level].push_back(
                this->multilevelcontext->getIndexofAboveCellContainingThisCell(level + 1, level, i));
          }
          // The above procedure might not be ordered and will create duplicates, get rid of them.
          sortAndEraseDuplicate(level);
        }
      }
    }

    void generateHierarchyBelowLevelExclusive(size_t coarsestLevel){

      // Do all level between this level and the finest
      for(size_t level=coarsestLevel + 1; level < this->multilevelcontext->getNumLevels() - 1; level++){

        auto currentLevelGrid = this->multilevelcontext->getGridForLevel(level);
        for(size_t i = 0; i < currentLevelGrid.size3; i++){
          size_t aboveindex = this->multilevelcontext->getIndexofAboveCellContainingThisCell(level, level - 1, i);
          if(this->isMasked(aboveindex, level - 1)){
            this->flaggedIdsAtEachLevel[level].push_back(i);
          }
        }
      }
    }

    void sortAndEraseDuplicate(size_t level){
      std::sort(this->flaggedIdsAtEachLevel[level].begin(), this->flaggedIdsAtEachLevel[level].end());
      this->flaggedIdsAtEachLevel[level].erase(std::unique(
          this->flaggedIdsAtEachLevel[level].begin(),
          this->flaggedIdsAtEachLevel[level].end()), this->flaggedIdsAtEachLevel[level].end());
    }

    bool isMasked(size_t id, size_t level){
      return std::binary_search(this->flaggedIdsAtEachLevel[level].begin(),
                                this->flaggedIdsAtEachLevel[level].end(), id);
    }
  };
}



#endif
