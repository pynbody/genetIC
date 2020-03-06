
#ifndef IC_EVALUATOR_HPP
#define IC_EVALUATOR_HPP

#include "../grid/grid.hpp"
#include "../grid/virtualgrid.hpp"
#include "field.hpp"
#include "../multilevelgrid/multilevelgrid.hpp"
#include <string>

namespace fields {
  /*! \class EvaluatorBase
      \brief The base class for field evaluators.
   *
   * Field evaluators exist to allow fields to be flexibly evaluated on grids which are not necessarily the same literal grid that
   * they are stored on. This allows the code to supersample or subsample a field onto a finer or coarser grid, for example.
   */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class EvaluatorBase : public std::enable_shared_from_this<EvaluatorBase<DataType, CoordinateType>> {
  public:
    //! \brief Returns the field evaluates at linear index i
    virtual DataType operator[](size_t i) const = 0;

    //! \brief Evaluates the field at the position 'at', using interpolation if necessary.
    virtual DataType operator()(const Coordinate<CoordinateType> &at) const = 0;

    //!\brief Returns true if index i corresponds to a point in the field.
    virtual bool contains(size_t i) const = 0;

    //! \brief Adds this field to the destination field.
    virtual void addTo(Field <DataType, CoordinateType> &destination) const {

      size_t size3 = destination.getGrid().size3;

      fields::cache::enableInterpolationCaches();

      // A note on the scheduling here: the most costly operations operations are encountered when there is
      // interpolation involved, at which point an thread-local LRU cache tries to reduce some of the workload. The hit
      // rate in the cache will be greatest if processors work on spatially localised areas. At present, such
      // strict localization is not implemented.

      /*
#pragma omp parallel for schedule(dynamic, destination.getGrid().size2)
      for (size_t ind_l = 0; ind_l < size3; ind_l++) {
        if (contains(ind_l))
          destination[ind_l] += (*this)[ind_l];
      }
       */

      destination.getGrid().parallelIterateOverCellsSpatiallyClustered([this, &destination](size_t ind_l) {
        if (contains(ind_l))
          destination[ind_l] += (*this)[ind_l];
      });


      fields::cache::disableInterpolationCaches();


    }
  };

  /*! \class DirectEvaluator
      \brief Evaluator that is appropriate when the grid and the field match perfectly.

  */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class DirectEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    const std::shared_ptr<const Field <DataType, CoordinateType>> field;

  public:
    DirectEvaluator(const Field <DataType, CoordinateType> &field) : field(field.shared_from_this()) {};

    //! \brief Direct evaluation of the field at a point where its value is stored.
    DataType operator[](size_t i) const override {
      return (*field)[i];
    }

    //! \brief Evaluate the field at a given point using interpolation
    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return field->evaluateInterpolated(at);
    }

    //! Check whether the specified point actually lies within this grid
    bool contains(size_t i) const override {
      return i < field->getGrid().size3;
    }
  };


  /*!   \class SuperSampleEvaluator
        \brief Evaluator that is appropriate when the grid is at a different resolution compared with the grid storage

        This evaluator generally always has to use interpolation to get the field at the specified points, because
        the data is not stored at the same resolution we are evaluating the field at.
  */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class SuperSampleEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::SuperSampleGrid<CoordinateType>;
    const std::shared_ptr<const grids::Grid<CoordinateType>> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying;

  public:
    SuperSampleEvaluator(const grids::VirtualGrid<CoordinateType> &grid,
                         std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying) :
      grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this())),
      underlying(underlying) {

    }

    //! \brief Interpolates to get the field at the centre of specified virtual cell.
    DataType operator[](size_t i) const override {
      auto centroid = grid->getCentroidFromIndex(i);
      return (*underlying)(centroid);
    }

    //! \brief Interpolates the underlying field at the given point
    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return (*underlying)(at);
    }

    bool contains(size_t i) const override {
      return grid->containsCell(i);
    }
  };

  /*!   \class SectionEvaluator
        \brief Evaluator that is appropriate when the grid is mapped using SectionOfGrid

        This corresponds to cases where the field is stored on a grid that expands beyond the confines of the grid on which we are currently evaluating.
  */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class SectionEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::SectionOfGrid<CoordinateType>;
    const std::shared_ptr<const grids::SectionOfGrid<CoordinateType>> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying;

  public:
    SectionEvaluator(const grids::VirtualGrid<CoordinateType> &grid,
                     std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying) :
      grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this())),
      underlying(underlying) {

    }

    //! Map to cell in underlying grid and evaluate directly there.
    DataType operator[](size_t i) const override {
      return (*underlying)[grid->mapIndexToUnderlying(i)];
    }

    //! \brief Evaluate the field at a given point using interpolation
    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return (*underlying)(at);
    }

    //! Check whether the supplied point actually lies in this SectionOfGrid
    bool contains(size_t i) const override {
      return grid->containsCell(i);
    }

  };

  /*!   \class SubSampleEvaluator
        \brief Evaluator used for sub-sampled grids.

        This class is used when we are accessing data at a lower resolution than it is stored.
        We have to use  coarse-graining to average over several stored field values to get the
        evaluated field.
  */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class SubSampleEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::SubSampleGrid<CoordinateType>;
    const std::shared_ptr<const grids::SubSampleGrid<CoordinateType>> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying;

  public:
    SubSampleEvaluator(const grids::VirtualGrid<CoordinateType> &grid,
                       std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlying) :
      grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this())),
      underlying(underlying) {

    }

    //! Average (coarse-grain) the points in the stored grid corresponding to the sub-sampled grid point
    virtual DataType operator[](size_t i) const override {
      DataType returnVal(0);
      CoordinateType localFactor3 = grid->forEachSubcell(i, [this, &returnVal](size_t local_id) {
        returnVal += (*(this->underlying))[local_id];
      });
      return returnVal / localFactor3;
    }

    //! \brief Evaluate the field at a given point using interpolation
    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      return (*underlying)(at);
    }

    bool contains(size_t i) const override {
      return grid->containsCell(i);
    }
  };

  /*!   \class ResolutionMatchingEvaluator
        \brief Evaluator used for resolution-matching grids.

        Resolution matching grids have a high-resolution window embedded in an interpolated low resolution box. Thus, this
        evaluator simply checks whether we are in the high resolution window or not, and passes the
        evaluation work to the relevant evaluator. It stores evaluators for both grids.
  */
  template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
  class ResolutionMatchingEvaluator : public EvaluatorBase<DataType, CoordinateType> {

  protected:
    using MyGridType = const grids::ResolutionMatchingGrid<CoordinateType>;
    const std::shared_ptr<MyGridType> grid;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlyingLoRes;
    const std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlyingHiRes;

  public:
    ResolutionMatchingEvaluator(const grids::VirtualGrid<CoordinateType> &grid,
                                std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlyingLoRes,
                                std::shared_ptr<const EvaluatorBase<DataType, CoordinateType>> underlyingHiRes) :

      grid(std::dynamic_pointer_cast<MyGridType>(grid.shared_from_this())),
      underlyingLoRes(underlyingLoRes), underlyingHiRes(underlyingHiRes) {

    }

    //! \brief Check whether i is in a high-resolution window and evaluate there if so; if not, directly evaluate it.
    virtual DataType operator[](size_t i) const override {
      auto coordinate = grid->getCoordinateFromIndex(i);
      if (grid->isInHiResWindow(coordinate)) {
        size_t mapped_index = grid->getIndexInHiResWindow(coordinate);
        return (*underlyingHiRes)[mapped_index];
      } else {
        return (*underlyingLoRes)[i];
      }
    }

    //! Interpolate, using the appropriate high resolution window if available for the given point.
    DataType operator()(const Coordinate<CoordinateType> &at) const override {
      if (grid->isInHiResWindow(at))
        return (*underlyingHiRes)(at);
      else
        return (*underlyingLoRes)(at);
    }

    bool contains(size_t i) const override {
      return grid->containsCell(i);
    }
  };


  //! \brief Return an object suitable for evaluating the specified field at coordinates on the specified grid
  template<typename DataType, typename CoordinateType>
  std::shared_ptr<EvaluatorBase<DataType, CoordinateType>> makeEvaluator(const Field <DataType, CoordinateType> &field,
                                                                         const grids::Grid<CoordinateType> &grid) {
                                                                         
    // Rather than duplicate the logic for a multi-level field (which is a more
    // general case), we create a dummy multi-level context and multi-level field
    // which in fact contain only a single level each, then call the multi-level
    // makeEvaluator.                     

    multilevelgrid::MultiLevelGrid<DataType> dummyContext;
    dummyContext.addLevel(const_cast<grids::Grid<CoordinateType> &>(field.getGrid()).shared_from_this());
    MultiLevelField<DataType> dummyMultiField(dummyContext,
                                              {const_cast<Field<DataType, CoordinateType> &>(field).shared_from_this()});
    return makeEvaluator(dummyMultiField, grid);

  };

  //! \brief Return an object suitable for evaluating the specified field at coordinates on the specified grid
  template<typename DataType, typename CoordinateType>
  std::shared_ptr<EvaluatorBase<DataType, CoordinateType>> makeEvaluator(const MultiLevelField <DataType> &field,
                                                                         const grids::Grid<CoordinateType> &grid) {

    // TODO: this routine, complete with use of RTTI, is really ugly and could do with being rethought.
    //
    // The overall strategy, which is probably sound, is first to check whether grid
    // is a base grid (rather than a virtual grid). If so, we locate the corresponding
    // single-level field within the multiple levels we have been provided with. (The
    // operation will throw an exception if there is no field defined for the grid
    // provided.) Otherwise, we create the appropriate adaptor evaluator and recurse
    // to find its underlying evaluators. Eventually this always bottoms out at a 
    // direct evaluator.

    auto &runtimeType = typeid(grid);

    if (runtimeType == typeid(grids::Grid<CoordinateType>)) {
      // Simplest case: the field is actually stored directly on this grid.
      return std::make_shared<DirectEvaluator<DataType, CoordinateType>>(field.getFieldForGrid(grid));
    } else {
      // Special case: ResolutionMatchingGrid points to TWO underlying grids
      if (runtimeType == typeid(grids::ResolutionMatchingGrid<CoordinateType>)) {
        const grids::ResolutionMatchingGrid<CoordinateType> &rmGrid =
          dynamic_cast<const grids::ResolutionMatchingGrid<CoordinateType> &>(grid);

        auto underlyingLoResGrid = rmGrid.getUnderlyingLoResInterpolated();
        auto underlyingHiResGrid = rmGrid.getUnderlyingHiRes();

        auto underlyingLoResEvaluator = makeEvaluator(field, *underlyingLoResGrid);
        auto underlyingHiResEvaluator = makeEvaluator(field, *underlyingHiResGrid);

        return std::make_shared<ResolutionMatchingEvaluator<DataType, CoordinateType>>
          (rmGrid, underlyingLoResEvaluator, underlyingHiResEvaluator);

      }

      // In all other cases, there is one underlying grid that we need the evaluator for, then we place an adaptor
      // around that.
      const grids::VirtualGrid<CoordinateType> &virtualGrid = dynamic_cast<const grids::VirtualGrid<CoordinateType> &>(grid);
      auto &underlyingGrid = *(virtualGrid.getUnderlying());
      auto underlyingEvaluator = makeEvaluator(field, underlyingGrid);

      if (runtimeType == typeid(grids::SectionOfGrid<CoordinateType>)) {
        return std::make_shared<SectionEvaluator<DataType, CoordinateType>>(virtualGrid, underlyingEvaluator);
      } else if (runtimeType == typeid(grids::SuperSampleGrid<CoordinateType>)) {
        return std::make_shared<SuperSampleEvaluator<DataType, CoordinateType>>(virtualGrid, underlyingEvaluator);
      } else if (runtimeType == typeid(grids::SubSampleGrid<CoordinateType>)) {
        return std::make_shared<SubSampleEvaluator<DataType, CoordinateType>>(virtualGrid, underlyingEvaluator);
      } else if (runtimeType == typeid(grids::OffsetGrid<CoordinateType>) ||
                 runtimeType == typeid(grids::MassScaledGrid<CoordinateType>) ||
                 runtimeType == typeid(grids::CenteredGrid<CoordinateType>) ||
                 runtimeType == typeid(grids::IndependentFlaggingGrid<CoordinateType>)) {
        return underlyingEvaluator;
      }
    }
    throw std::runtime_error(std::string("Don't know how to evaluate field on grid of type ") + runtimeType.name());

  };


}

#endif
