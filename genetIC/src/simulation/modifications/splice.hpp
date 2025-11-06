#ifndef IC_SPLICE_HPP
#define IC_SPLICE_HPP

#include <complex>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/numerics/cg.hpp>
#include <src/tools/numerics/minres.hpp>

namespace modifications {
  template<typename T>
  fields::Field<char, T> generateMaskFromFlags(const grids::Grid<T> &grid) {
    std::vector<size_t> flags;
    grid.getFlaggedCells(flags);

    fields::Field<char, T> mask(const_cast<grids::Grid<T> &>(grid), false); // fills with zeros / falses
    for(auto f: flags) {
      mask[f] = true;
    }
    return mask;
  }

  template<typename T>
  fields::Field<char, T> generateMaskComplementFromFlags(const grids::Grid<T> &grid) {
    std::vector<size_t> flags;
    grid.getFlaggedCells(flags);

    fields::Field<char, T> mask(const_cast<grids::Grid<T> &>(grid), false); // fills with zeros / falses
    for(size_t i=0; i<mask.getDataVector().size(); ++i) {
      mask[i] = true;
    }
    for(auto f: flags) {
      mask[f] = false;
    }
    return mask;
  }

  // Compute the operator that applies T^+ ... T *at all levels*
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> T_op_T (
      const fields::OutputField<DataType> &inputs,
      const std::vector<fields::Field<DataType, T>> &covs,
      const std::function<void(const int, fields::Field<DataType, T> &)> op
  ) {
    fields::OutputField<DataType> outputs(inputs.getContext(), inputs.getTransferType());
    outputs.getFieldForLevel(0); // trigger allocation

    fields::OutputField<DataType> inputs_delta(inputs);
    auto filters = inputs.getFilters();
    int Nlevel = inputs.getNumLevels();

    const auto& multiLevelContext = outputs.getContext();

    assert (inputs.isRealOnAllLevels());

    // Apply C^0.5 to the input fields
    for (size_t level=0; level<Nlevel; ++level) {
      auto& field = inputs_delta.getFieldForLevel(level);
      field.toFourier();
      field.applyTransferFunction(covs[level], 0.5);
    }

    // Filter all levels (but the last) in their windows
    // Notice we are doing window then filter
    for (size_t level=0; level<Nlevel; ++level) {
      auto& field = inputs_delta.getFieldForLevel(level);
      const auto &f = filters.getFilterForLevel(level);
      if (level < Nlevel - 1) {
        auto window = multiLevelContext.getGridForLevel(level+1).getWindow();
        field.applyFilterInWindow(f, window, true);
      } else {
        field.toFourier();
        field.applyFilter(f);
      }
    }

    // --------------------------------------------------------------------
    // Operator applies to the field + contributions from all finer levels
    for (size_t level=0; level<Nlevel; ++level) {
      auto out = inputs_delta.getFieldForLevel(level).copy();

      // Add contributions from all coarser levels
      for (size_t source_level = 0; source_level < level; ++source_level) {
        auto source_field = inputs_delta.getFieldForLevel(source_level).copy();
        T pixel_volume_ratio = multiLevelContext.getWeightForLevel(level) /
                               multiLevelContext.getWeightForLevel(source_level);
        
        source_field->toReal();
        *source_field *= sqrt(pixel_volume_ratio);

        out->addFieldFromDifferentGrid(*source_field);
      }

      // Apply operator
      op(level, *out);

      // Add contributions from all finer levels
      for (size_t source_level = level + 1; source_level < Nlevel; ++source_level) {
        auto source_field = inputs_delta.getFieldForLevel(source_level).copy();
        T pixel_volume_ratio = multiLevelContext.getWeightForLevel(level) /
                               multiLevelContext.getWeightForLevel(source_level);

        source_field->toFourier();
        *source_field *= sqrt(pixel_volume_ratio);
        op(source_level, *source_field);

        out->addFieldFromDifferentGrid(*source_field);
      }

      outputs.getFieldForLevel(level) = std::move(*out);
    }

    // Filter all levels (but the last) in their windows
    // Notice we are doing filter then window (the opposite order to above)
    for (size_t level=0; level<Nlevel; ++level) {
      auto& field = outputs.getFieldForLevel(level);
      const auto &f = filters.getFilterForLevel(level);
      if (level < Nlevel - 1) {
        auto window = multiLevelContext.getGridForLevel(level+1).getWindow();
        field.applyFilterInWindow(f, window, false);
      } else {
        field.toFourier();
        field.applyFilter(f);
      }
    }
    
    // Apply C^0.5 to the input fields
    for (size_t level=0; level<Nlevel; ++level) {
      auto& field = outputs.getFieldForLevel(level);
      field.toFourier();
      field.applyTransferFunction(covs[level], 0.5);
    }

    outputs.toReal();
    return outputs;
  };

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> Mbar_Cm1_Mbar(
    const fields::OutputField<DataType> & inputs,
    const auto& covs, const auto& masks, const auto& masksCompl, size_t Nlevel
  ) {
    auto outputs = T_op_T<DataType, T>(
      inputs,
      covs,
      [&](const int level, fields::Field<DataType,T> & input) -> void
    {
      input.toReal();
      input *= masksCompl[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= masksCompl[level];
    });
    return outputs;
  }

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> Mbar_Cm1_M(
    const fields::OutputField<DataType> & inputs,
    const auto& covs, const auto& masks, const auto& masksCompl, size_t Nlevel
  ) {
    auto outputs = T_op_T<DataType, T>(
      inputs, covs,
      [&](const int level, fields::Field<DataType,T> & input) -> void
    {
      input.toReal();
      input *= masks[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= masksCompl[level];
    });
    return outputs;
  };

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> splice(fields::OutputField<DataType> & a,
                                       fields::OutputField<DataType> & b) {

      assert (a.getTransferType() == particle::species::whitenoise);
      assert (b.getTransferType() == particle::species::whitenoise);
      assert (a.isFourierOnAllLevels());
      assert (b.isFourierOnAllLevels());

      // To understand the implementation below, first read Appendix A of Cadiou et al (2021),
      // and/or look at the 1D toy implementation (in tools/toy_implementation/gene_splicing.ipynb) which
      // contains a similar derivation and near-identical implementation.

      std::vector<fields::Field<DataType,T>> covs;
      std::vector<fields::Field<char,T>> masks;
      std::vector<fields::Field<char,T>> masksCompl;

      int Nlevel = a.getNumLevels();
      for(size_t level=0; level<Nlevel; ++level) {
        auto ctxt = a.getContext();
        fields::Field<DataType,T> cov(*ctxt.getCovariance(level, particle::species::all));
        cov.setFourierCoefficient(0, 0, 0, 1);
        covs.push_back(cov);

        masks.push_back(generateMaskFromFlags(ctxt.getGridForLevel(level)));
        masksCompl.push_back(generateMaskComplementFromFlags(ctxt.getGridForLevel(level)));
      }

      a.toReal();
      b.toReal();

      fields::OutputField<T> delta(b);
      delta -= a;

      fields::OutputField<T> z = Mbar_Cm1_M(delta, covs, masks, masksCompl, Nlevel);

      auto A = [&](const fields::OutputField<T> & inputs) {
        return Mbar_Cm1_Mbar(inputs, covs, masks, masksCompl, Nlevel);
      };


      // fields::OutputField<DataType> alpha = tools::numerics::conjugateGradient<DataType>(A, z);
      fields::OutputField<DataType> alpha = tools::numerics::minres<DataType>(A, z);

      // alpha.toFourier();
      // alpha.applyTransferFunction(preconditioner, 0.5);
      // alpha.toReal();

      // fields::Field<DataType,T> bInDeltaBasis(b);
      // bInDeltaBasis.toFourier();
      // bInDeltaBasis.applyTransferFunction(preconditioner, 0.5);
      // bInDeltaBasis.toReal();

      // alpha*=maskCompl;
      // alpha+=bInDeltaBasis;

      // delta_diff*=mask;
      // alpha-=delta_diff;

      // assert(!alpha.isFourier());
      // alpha.toFourier();
      // alpha.applyTransferFunction(preconditioner, -0.5);
      // alpha.toReal();

      return alpha;
  }
}

#endif //IC_SPLICE_HPP
