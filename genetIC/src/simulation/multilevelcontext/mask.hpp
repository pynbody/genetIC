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

  protected:
    MultiLevelContextInformation<DataType>* multilevelcontext;
    std::vector<std::vector<size_t>> flaggedIdsAtEachLevel;

    virtual void generateFlagsHierarchy() = 0;

  };

  template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
  class RamsesMask: public AbstractBaseMask<DataType, T> {

  public:
    explicit RamsesMask(MultiLevelContextInformation<DataType>* multilevelcontext_):
        AbstractBaseMask<DataType, T> (multilevelcontext_){
      generateFlagsHierarchy();
      ensureFlaggedVolumeIsContinuous();
    };

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

    void ensureFlaggedVolumeIsContinuous() override{
      //TODO Add Lagrangian volume definition if it bothers you
    }

    void recalculateWithNewContext(MultiLevelContextInformation<DataType>* multilevelcontext_){
      this->multilevelcontext = multilevelcontext_;
      this->flaggedIdsAtEachLevel.clear();
      this->flaggedIdsAtEachLevel = std::vector<std::vector<size_t>>(this->multilevelcontext->getNumLevels());
      generateFlagsHierarchy();
      ensureFlaggedVolumeIsContinuous();
    }

  protected:
    void generateFlagsHierarchy() override {
      // From deepest flagged level, generate vector of flagged cells at each level
      size_t deepestFlaggedLevel = this->multilevelcontext->deepestLevelwithFlaggedCells();

      for (int level = deepestFlaggedLevel; level >= 0; level--) {

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
          std::sort(this->flaggedIdsAtEachLevel[level].begin(), this->flaggedIdsAtEachLevel[level].end());
          this->flaggedIdsAtEachLevel[level].erase(std::unique(
              this->flaggedIdsAtEachLevel[level].begin(),
              this->flaggedIdsAtEachLevel[level].end()), this->flaggedIdsAtEachLevel[level].end());
        }
      }
    }
  };
}



#endif
