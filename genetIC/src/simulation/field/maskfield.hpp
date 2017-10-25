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

  };

  //! This implementation of mask is targeted for RAMSES refinement masks.
  template<typename DataType>
  class RAMSESMaskField : public AbstractMaskField<DataType> {

  public:
    explicit RAMSESMaskField(multilevelcontext::MultiLevelContextInformation<DataType> &multiLevelContext)
        : AbstractMaskField<DataType>(multiLevelContext){

    }

    void calculateMasksAllLevels() override {
      size_t finest_level = this->multiLevelContext->getNumLevels() - 1;

      if (finest_level == 0) {
        std::cerr << "GRAFIC: Not generating a mask for a uniform resolution grid" << std::endl;
      } else {
        calculateMaskFinestLevel();

        for (size_t level = 0; level < finest_level - 1; ++level) {
          auto current_level_grid = this->multiLevelContext->getGridForLevel(level);
          auto next_level_grid = this->multiLevelContext->getGridForLevel(level + 1);
          refineNextZoomGrid(current_level_grid, next_level_grid, level);
        }
      }
    }

  protected:
    using T = typename AbstractMaskField<DataType>::T;

    void calculateMaskFinestLevel() override {
      size_t finest_level = this->multiLevelContext->getNumLevels() - 1;
      size_t deepest_flagged_level = finest_level;
      // There might be several virtual grids above with no flagged cells to refer to
      //TODO Unfinished, need to catch runtime err if there are no flagged cells at all
      try {
        deepest_flagged_level = this->multiLevelContext->deepestLevelwithFlaggedCells();
      } catch (std::runtime_error& err) {

      }

      if (deepest_flagged_level == finest_level) {
        throw std::runtime_error("There are some flagged cells on the finest level. This is weird...");
      }

      std::vector<size_t> deepest_flagged_cells;
      this->multiLevelContext->getGridForLevel(deepest_flagged_level).getFlaggedCells(deepest_flagged_cells);

      refineFlaggedCellswrtoAbove(this->multiLevelContext->getGridForLevel(finest_level),
                                  this->multiLevelContext->getGridForLevel(deepest_flagged_level),
                                  finest_level);

    };

  private:
    void refineNextZoomGrid(grids::Grid<DataType>& current_level_grid,
                                 grids::Grid<DataType>& other_level_grid,
                                 size_t level){
#pragma omp parallel for
      for (size_t i = 0; i < current_level_grid.size3; i++) {
        Coordinate<T> cell_coord(current_level_grid.getCellCentroid(i));
        if (other_level_grid.containsPoint(cell_coord)) {
          this->getFieldForLevel(level).getDataVector()[i] = 1;
        }
      }
    }

    void refineFlaggedCellsOnSameLevel(std::vector<size_t>& flaggedcellsarray,
                              size_t level){
#pragma omp parallel for
      for (size_t i=0; i< flaggedcellsarray.size(); i++) {
        this->getFieldForLevel(level).getDataVector()[flaggedcellsarray[i]] = 1;
      }
    }


    void refineFlaggedCellswrtoAbove(grids::Grid<DataType>& current_level_grid,
                                     grids::Grid<DataType>& other_level_grid,
                                     size_t level){

      std::vector<size_t> flags;
      other_level_grid.getFlaggedCells(flags);

#pragma omp parallel for
      for (size_t i = 0; i < current_level_grid.size3; i++) {
        Coordinate<T> cell_coord(current_level_grid.getCellCentroid(i));
         size_t id_on_other = other_level_grid.getCellContainingPoint(cell_coord);

        // If the current cell is flagged on the above level, mark it valid for refinement
        if((std::binary_search(flags.begin(), flags.end(), id_on_other))){
          this->getFieldForLevel(level).getDataVector()[i] = 1;
        }
      }
    }

    void refineEntireGrid(grids::Grid<DataType>& current_level_grid, size_t level){
#pragma omp parallel for
      for (size_t i = 0; i < current_level_grid.size3; i++) {
        this->getFieldForLevel(level).getDataVector()[i] = 1;
      }
    }
};





}
#endif
