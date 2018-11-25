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

  //! Extend the concept of mask for grafic outputs. The masks must be carried through the entire grafic hierarchy, including virtual intermediate grids.
  //! Useful for generating ic_refmap/ic_pvar files for RAMSES for example
  template<typename DataType, typename T=tools::datatypes::strip_complex <DataType>>
  class GraficMask: public AbstractBaseMask<DataType, T> {

  private:
    const std::vector<std::vector<size_t>> &inputzoomParticlesAsMask;
    int deepestLevelWithMaskedCells=-1;

  public:
    explicit GraficMask(MultiLevelContextInformation<DataType>* multilevelcontext_, std::vector<std::vector<size_t>> &input_mask):
        AbstractBaseMask<DataType, T> (multilevelcontext_), inputzoomParticlesAsMask(input_mask){
      assert(inputzoomParticlesAsMask.size() <= this->flaggedIdsAtEachLevel.size());
    };

    /*!
     * @return 1.0 if cell with index is in the mask in this level, 0.0 else.
     */
    T const isInMask(size_t level, size_t cellindex) override{

      if(this->flaggedIdsAtEachLevel[level].size() == 0){
        // No flagged cells at all on this level, refine everywhere
        return T(1.0);
      } else if(isMasked(cellindex, level)){
        // Cell is among flagged cells, refine it
        return T(1.0);
      }
      return T(0.0);

    }

    void calculateMask() override{
      identifyLevelsOfInputMask();
      generateFlagsHierarchy();
      ensureFlaggedVolumeIsContinuous();
    }

    void ensureFlaggedVolumeIsContinuous() override{
      //TODO Add Lagrangian volume definition if it bothers you
    }


  protected:
    //! Information in the input mask needs to be matched to a full Grafic hierarchy, potentially with virtual intermediate levels.
    void identifyLevelsOfInputMask() {
      size_t input_level=0;

      for (size_t level = 0; level < this->multilevelcontext->getNumLevels() - 1; level++) {
        for (size_t i = 0; i < this->multilevelcontext->getGridForLevel(level).size3; i++) {

          // Terminate calculation if all possible input masks have been found
          if(input_level >= this->inputzoomParticlesAsMask.size()){
            return;
          }

          // As soon as one cell id matches, we have identified the level, so assign it.
          if (this->isPartofInputMask(i, input_level)) {

            assert(input_level < this->inputzoomParticlesAsMask.size());
            assert(int(level)>= this->deepestLevelWithMaskedCells);

            this->flaggedIdsAtEachLevel[level] = this->inputzoomParticlesAsMask[input_level];
            this->deepestLevelWithMaskedCells = int(level);
            input_level++;
            break;
          }
        }
      }
    }

    //! Calculates the flagged cells on all levels
    /*!
     * For each level, imports the flagged cells of the level if there are any. If not, e.g.
     * for virtual grids, downgrades the higher resolution mask on the virtual grid.
     */
    void generateFlagsHierarchy() override {

      if(this->deepestLevelWithMaskedCells<0){
        std::cerr << "WARNING No zoom regions were ever opened. Grafic mask will not be generated in this case" <<std::endl;
        return;
      }

      generateHierarchyAboveLevelInclusive(this->deepestLevelWithMaskedCells);
      generateHierarchyBelowLevelExclusive(this->deepestLevelWithMaskedCells);
    }

    void generateHierarchyAboveLevelInclusive(int deepestLevel){
      for (int level = deepestLevel; level >= 0; level--) {

        if (this->flaggedIdsAtEachLevel[level].size() == 0) {
        // Generate flags on intermediate levels that will not have some
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
    }

    bool isMasked(size_t id, size_t level){
      return std::binary_search(this->flaggedIdsAtEachLevel[level].begin(),
                                this->flaggedIdsAtEachLevel[level].end(), id);
    }

    bool isPartofInputMask(size_t id, size_t level){
      return std::binary_search(this->inputzoomParticlesAsMask[level].begin(),
                                this->inputzoomParticlesAsMask[level].end(), id);
    }

  public:
    //! Generate a full multilevel field storing the mask information.
    /*! Increases memory requirement for the code but useful for debugging.
     */
    std::shared_ptr<fields::MultiLevelField<DataType>> convertToField(){
      std::vector<std::shared_ptr<fields::Field<DataType, T>>> fields;

      // Field full of zeros
      for (size_t level = 0; level < this->multilevelcontext->getNumLevels(); ++level) {
        fields.push_back(std::shared_ptr<fields::Field<DataType, T>>(
            new fields::Field<DataType, T>(this->multilevelcontext->getGridForLevel(level))));
      }

      auto maskfield = std::make_shared<fields::MultiLevelField<DataType>>(*(this->multilevelcontext), fields);

      // Field with mask information
      for (size_t level = 0; level < this->multilevelcontext->getNumLevels(); ++level) {
        for (size_t i_z = 0; i_z < this->multilevelcontext->getGridForLevel(level).size; ++i_z) {
          for (size_t i_y = 0; i_y < this->multilevelcontext->getGridForLevel(level).size; ++i_y) {
            for (size_t i_x = 0; i_x < this->multilevelcontext->getGridForLevel(level).size; ++i_x) {

              // These two indices can be different if some virtual grid are used in the context, e.g. centered.
              // In all other cases, they will be equal.
              size_t i = size_t(i_x * this->multilevelcontext->getGridForLevel(level).size + i_y)
                         * this->multilevelcontext->getGridForLevel(level).size + i_z;
              size_t virtual_i = this->multilevelcontext->getGridForLevel(level).getIndexFromCoordinateNoWrap(i_x, i_y, i_z);

              maskfield->getFieldForLevel(level).getDataVector()[i] = isInMask(level, virtual_i);

            }
          }
        }
      }
      return maskfield;
    }
  };
}



#endif
