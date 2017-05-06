
#ifndef IC_EVALUATOR_HPP
#define IC_EVALUATOR_HPP

#include "../grid/grid.hpp"
#include "../grid/virtualgrid.hpp"
#include "field.hpp"
#include <string>

namespace fields {
  /*! The base class for field evaluators.
   *
   * Field evaluators exist to allow fields to be evaluated on grids which are not the same literal grid that
   * they are stored on. This is useful for example when supersampling or subsampling.
   */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class EvaluatorBase : public std::enable_shared_from_this<EvaluatorBase<DataType, CoordinateType>> {
  public:
    virtual DataType operator[](size_t i) const =0;
    virtual DataType operator()(const Coordinate<CoordinateType> &at) const = 0;
    virtual bool contains(size_t i) const =0;

    void addTo(Field<DataType, CoordinateType> &destination) {

      size_t size3 = destination.getGrid().size3;

// #pragma omp parallel for schedule(static)
      for (size_t ind_l = 0; ind_l < size3; ind_l++) {
        if(contains(ind_l))
          destination[ind_l] += (*this)[ind_l];
      }
    }
  };

  /*! Evaluator that is appropriate when the grid and the field match perfectly.
   *
   */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class DirectEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    const std::shared_ptr<const Field<DataType, CoordinateType>> field;

  public:
    DirectEvaluator(const Field<DataType, CoordinateType> & field) : field(field.shared_from_this()) {};

    DataType operator[](size_t i) const override {
      return (*field)[i];
    }

    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return field->evaluateInterpolated(at);
    }

    bool contains(size_t i) const override {
      return i<field->getGrid().size3;
    }
  };


  //! Evaluator that is appropriate when the grid is at a different resolution compared with the grid storage
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class SuperSampleEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::SuperSampleGrid<CoordinateType>;
    const std::shared_ptr<const grids::Grid<CoordinateType>> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying;

  public:
    SuperSampleEvaluator(std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying,
    const grids::VirtualGrid<CoordinateType> &grid) :
    underlying(underlying), grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this()))
    {

    }

    DataType operator[](size_t i) const override {
      auto centroid = grid->getCellCentroid(i);
      return (*underlying)(centroid);
    }

    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return (*underlying)(at);
    }

    bool contains(size_t i) const override {
      return i<grid->size3;
    }
  };

  //! Evaluator that is appropriate when the grid is mapped using SectionOfGrid
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class SectionEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::SectionOfGrid<CoordinateType>;
    const std::shared_ptr<const grids::SectionOfGrid<CoordinateType>> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying;

  public:
    SectionEvaluator(std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying,
                     const grids::VirtualGrid<CoordinateType> &grid) :
        underlying(underlying), grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this()))
    {

    }

    DataType operator[](size_t i) const override {
      return (*underlying)[grid->mapIndexToUnderlying(i)];
    }

    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return (*underlying)(at);
    }

    bool contains(size_t i) const override {
      return grid->containsCell(i);
    }

  };


  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class SubSampleEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::SubSampleGrid<CoordinateType>;
    const std::shared_ptr<const grids::SubSampleGrid<CoordinateType>> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying;

  public:
    SubSampleEvaluator(std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying,
    const grids::VirtualGrid<CoordinateType> &grid) :
    underlying(underlying), grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this()))
    {

    }

    virtual DataType operator[](size_t i) const override {
      DataType returnVal(0);
      int localFactor3 = grid->forEachSubcell(i, [this, &returnVal](size_t local_id) {
        returnVal += (*(this->underlying))[local_id];
      });
      return returnVal / localFactor3;
    }

    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return (*underlying)(at);
    }

    bool contains(size_t i) const override {
      return i<grid->size3;
    }
  };



  template<typename DataType, typename CoordinateType>
  std::shared_ptr<EvaluatorBase<DataType, CoordinateType>> makeEvaluator(const Field<DataType, CoordinateType> &field,
                                               const grids::Grid<CoordinateType> &grid) {
    if(!grid.pointsToGrid(&field.getGrid())) {
      throw std::runtime_error("Cannot evaluate this Field on the specified Grid because the grids do not have a well-defined relationship to each other.");
    }

    auto &runtimeType = typeid(grid);

    if(runtimeType==typeid(grids::Grid<CoordinateType>)) {
      return std::make_shared<DirectEvaluator<DataType, CoordinateType>>(field);
    } else {
      const grids::VirtualGrid<CoordinateType> &virtualGrid = dynamic_cast<const grids::VirtualGrid<CoordinateType> &>(grid);
      auto & underlyingGrid = *(virtualGrid.getUnderlying());
      auto underlyingEvaluator = makeEvaluator(field, underlyingGrid);

      if(runtimeType==typeid(grids::SectionOfGrid<CoordinateType>)) {
        return std::make_shared<SectionEvaluator<DataType, CoordinateType>>(underlyingEvaluator, virtualGrid);
      } else if(runtimeType==typeid(grids::SuperSampleGrid<CoordinateType>)) {
        return std::make_shared<SuperSampleEvaluator<DataType, CoordinateType>>(underlyingEvaluator, virtualGrid);
      } else if(runtimeType==typeid(grids::SubSampleGrid<CoordinateType>)) {
        return std::make_shared<SubSampleEvaluator<DataType, CoordinateType>>(underlyingEvaluator, virtualGrid);
      } else if(runtimeType==typeid(grids::OffsetGrid<CoordinateType>) ||
          runtimeType==typeid(grids::MassScaledGrid<CoordinateType>)) {
        return underlyingEvaluator;
      }
    }
    throw std::runtime_error(std::string("Don't know how to evaluate field on grid of type ")+runtimeType.name());
  };



}

#endif //IC_EVALUATOR_HPP
