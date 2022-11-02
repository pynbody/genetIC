#ifndef IC_SPLICE_HPP
#define IC_SPLICE_HPP

#include <complex>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/numerics/cg.hpp>

namespace modifications {

  enum TOPERATOR {
    TPLUS, TDAGGER
  };

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

  //! Return the field f which satisfies f = a in flagged region while minimising (f-b).C^-1.(f-b) elsewhere
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::Field<DataType,T> spliceOneLevel(fields::Field<DataType,T> & a,
                                           fields::Field<DataType,T> & b,
                                           const fields::Field<DataType,T> & cov) {

      // To understand the implementation below, first read Appendix A of Cadiou et al (2021),
      // and/or look at the 1D toy implementation (in tools/toy_implementation/gene_splicing.ipynb) which
      // contains a similar derivation and near-identical implementation.

      assert (&a.getGrid() == &b.getGrid());
      assert (&a.getGrid() == &cov.getGrid());
      auto mask = generateMaskFromFlags(a.getGrid());
      auto maskCompl = generateMaskComplementFromFlags(a.getGrid());

      a.toFourier();
      b.toFourier();

      // The preconditioner should be almost equal to the covariance.
      // We however set the fundamental of the power spectrum to a non-null value,
      // otherwise, the mean value in the spliced region is unconstrained.
      fields::Field<DataType,T> preconditioner(cov);
      preconditioner.setFourierCoefficient(0, 0, 0, 1);

      fields::Field<DataType,T> delta_diff = b-a;
      delta_diff.applyTransferFunction(preconditioner, 0.5);
      delta_diff.toReal();

      fields::Field<DataType,T> z(delta_diff);
      z*=mask;
      z.toFourier();
      z.applyTransferFunction(preconditioner, -1.0);
      z.toReal();
      z*=maskCompl;
      z.toFourier();
      z.applyTransferFunction(preconditioner, 0.5);
      z.toReal();

      auto X = [&preconditioner, &maskCompl](const fields::Field<DataType,T> & input) -> fields::Field<DataType,T>
      {
        fields::Field<DataType,T> v(input);
        assert (!input.isFourier());
        assert (!v.isFourier());
        v.toFourier();
        v.applyTransferFunction(preconditioner, 0.5);
        v.toReal();
        v*=maskCompl;
        v.toFourier();
        v.applyTransferFunction(preconditioner, -1.0);
        v.toReal();
        v*=maskCompl;
        v.toFourier();
        v.applyTransferFunction(preconditioner, 0.5);
        v.toReal();
        return v;
      };

      fields::Field<DataType,T> alpha = tools::numerics::conjugateGradient<DataType>(X, z);

      alpha.toFourier();
      alpha.applyTransferFunction(preconditioner, 0.5);
      alpha.toReal();

      fields::Field<DataType,T> bInDeltaBasis(b);
      bInDeltaBasis.toFourier();
      bInDeltaBasis.applyTransferFunction(preconditioner, 0.5);
      bInDeltaBasis.toReal();

      alpha*=maskCompl;
      alpha+=bInDeltaBasis;

      delta_diff*=mask;
      alpha-=delta_diff;

      assert(!alpha.isFourier());
      alpha.toFourier();
      alpha.applyTransferFunction(preconditioner, -0.5);
      alpha.toReal();

      return alpha;
  }

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> spliceMultiLevel(
    const fields::OutputField<DataType> &noiseB,
    const fields::OutputField<DataType> &noiseA
  ) {
    auto multiLevelContext = noiseA.getContext();

    // Prefetch covariances and masks
    int Nlevel = noiseA.getNumLevels();
    std::vector<fields::Field<char, T>> masksCompl;
    std::vector<fields::Field<char, T>> masks;
    std::vector<fields::Field<DataType, T>> covs;

    auto filters = noiseA.getFilters();

    assert (noiseA.getTransferType() == particle::species::whitenoise);
    assert (noiseB.getTransferType() == particle::species::whitenoise);


    for(size_t level=0; level<Nlevel; ++level) {
      masksCompl.push_back(generateMaskComplementFromFlags(multiLevelContext.getGridForLevel(level)));
      masks.push_back(generateMaskFromFlags(multiLevelContext.getGridForLevel(level)));

      fields::Field<DataType, T> cov(*multiLevelContext.getCovariance(level, particle::species::all));
      cov.setFourierCoefficient(0, 0, 0, 1);
      covs.push_back(cov);
    }

    // Compute the operator that applies T^+ ... T *at all levels*
    auto T_op_T = [&covs, &filters, &multiLevelContext, Nlevel](
        const fields::OutputField<DataType> &inputs,
        const std::function<void(const int, fields::Field<DataType, T> &)> op,
        const TOPERATOR Top
    ) -> fields::OutputField<DataType> {
      fields::OutputField<DataType> outputs(inputs.getContext(), inputs.getTransferType());
      fields::OutputField<DataType> inputCopy(inputs);
      const double power_out = (Top == TDAGGER) ? 0.5 : -0.5;

      assert (inputs.isRealOnAllLevels());

      // Compute operator using each possible pair
      for (auto level=0; level<Nlevel; ++level) {
        auto & out = outputs.getFieldForLevel(level);
        out.toReal();

        for (auto level_other=std::max(0, level-1); level_other<std::min(Nlevel, level+2); ++level_other) {
          T pixel_volume_ratio = multiLevelContext.getWeightForLevel(level) /
                                 multiLevelContext.getWeightForLevel(level_other);

          fields::Field<DataType, T> tmp(out.getGrid(), false);

          if (level_other > level) {        // Contribution from finer level
            auto tmp_from_other = inputCopy.getFieldForLevel(level_other).copy();
            tmp_from_other->toReal();
            tmp.toReal();
            tmp.addFieldFromDifferentGrid(*tmp_from_other);
            tmp.toFourier();
            tmp.applyFilter(filters.getFilterForLevel(level));
            tmp.applyTransferFunction(covs[level], 0.5);
            tmp.toReal();
            op(level, tmp);
            tmp.toFourier();
            tmp.applyTransferFunction(covs[level], power_out);
            tmp.applyFilter(filters.getFilterForLevel(level));
            tmp.toReal();
          } else if (level_other < level) { // Contribution from coarser level
            auto tmp_from_other = inputCopy.getFieldForLevel(level_other).copy();
            tmp_from_other->toFourier();
            tmp_from_other->applyFilter(filters.getFilterForLevel(level_other));
            tmp_from_other->applyTransferFunction(covs[level_other], 0.5);
            tmp_from_other->toReal();
            op(level_other, *tmp_from_other);
            tmp_from_other->toFourier();
            tmp_from_other->applyTransferFunction(covs[level_other], power_out);
            tmp_from_other->applyFilter(filters.getFilterForLevel(level_other));
            tmp_from_other->toReal();
            tmp.toReal();
            tmp.addFieldFromDifferentGrid(*tmp_from_other);
          } else {
            tmp += inputCopy.getFieldForLevel(level);
            tmp.toFourier();
            tmp.applyFilter(filters.getFilterForLevel(level));
            tmp.applyTransferFunction(covs[level], 0.5);
            tmp.toReal();
            op(level, tmp);
            tmp.toFourier();
            tmp.applyTransferFunction(covs[level], power_out);
            tmp.applyFilter(filters.getFilterForLevel(level));
            tmp.toReal();
          }
          out.addScaled(tmp, std::sqrt(pixel_volume_ratio));
        }
      }

      return outputs;
    };

    auto Tplus_op_T = [&T_op_T](const fields::OutputField<DataType> &inputs,
                                const std::function<void(const int, fields::Field<DataType, T> &)> op
    ) -> fields::OutputField<DataType> {
      return T_op_T(inputs, op, TPLUS);
    };

    auto Tdagger_op_T = [&T_op_T](const fields::OutputField<DataType> &inputs,
                                  const std::function<void(const int, fields::Field<DataType, T> &)> op
    ) -> fields::OutputField<DataType> {
      return T_op_T(inputs, op, TDAGGER);
    };

    // Compute Mbar C^-1 Mbar for a single level
    auto Mbar_Cm1_Mbar_op = [&covs, &multiLevelContext, &masksCompl](const int level, fields::Field<DataType,T> & input) {
      input.toReal();
      input *= masksCompl[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= masksCompl[level];
    };

    // Compute Mbar C^-1 M for a single field
    auto Mbar_Cm1_M_op = [&covs, &multiLevelContext, &masksCompl, &masks](const int level, fields::Field<DataType,T> & input) -> fields::Field<DataType,T> {
      input.toReal();
      input *= masks[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= masksCompl[level];
      return input;
    };

    auto preconditioner = [&covs, &masksCompl, &filters, Nlevel, &Tdagger_op_T, &Tplus_op_T](fields::OutputField<DataType> & inputs) {
      // const bool precondition = false;

      // if (precondition) {
      //   for (auto level=0; level<Nlevel; ++level) {
      //     auto & field = inputs.getFieldForLevel(level);
      //     field.toFourier();
      //     field.applyTransferFunction(covs[level], -0.5);
      //     field.applyFilter(filters.getFilterForLevel(level));
      //     field.toReal();
      //     field *= masksCompl[level];
      //     field.toFourier();
      //     field.applyTransferFunction(covs[level], 1.0);
      //     field.toReal();
      //     field *= masksCompl[level];
      //     field.toFourier();
      //     field.applyTransferFunction(covs[level], -0.5);
      //     field.applyFilter(filters.getFilterForLevel(level));
      //     field.toReal();
      //   }
      // }
    };

    // Full splicing operator
    // Q := T^+ Mbar C^-1 Mbar T
    auto Q = [&Tdagger_op_T, &Tplus_op_T, &Mbar_Cm1_Mbar_op, &preconditioner](
      const fields::OutputField<DataType> & inputs
    ) -> fields::OutputField<DataType> {
      auto outputs = inputs;
      preconditioner(outputs);
      outputs = Tdagger_op_T(outputs, Mbar_Cm1_Mbar_op);
      preconditioner(outputs);
      outputs.toReal();
      return outputs;
    };

    // Compute the rhs T^+ Mbar C^-1 M T (b - a)
    auto rhs = noiseB;
    rhs -= noiseA;
    // rhs.toFourier();
    // rhs.applyPowerSpectrumFor(particle::species::all);
    // rhs.combineGrids();
    // rhs.applyPowerSpectrumFor(particle::species::whitenoise);
    rhs.toReal();
    rhs = Tdagger_op_T(rhs, Mbar_Cm1_M_op);
    preconditioner(rhs);
    rhs.toReal();

    // Solve the linear system [T^t Mbar C^-1 Mbar T] n_alpha = rhs
    // that approximates [T^+ Mbar C^-1 Mbar T] and is symmetric
    auto n_alpha = tools::numerics::minres<DataType>(Q, rhs);
    // auto n_alpha = tools::numerics::bicgstab<DataType>(Q, rhs);
    preconditioner(n_alpha);


    // T^+ M T (n_a-n_b)
    auto Ma_minus_b = noiseA;
    Ma_minus_b -= noiseB;
    Ma_minus_b.toFourier();
    Ma_minus_b.applyPowerSpectrumFor(particle::species::all);
    Ma_minus_b.combineGrids();
    Ma_minus_b.toReal();
    for (auto level = 0; level<Nlevel; ++level)
      Ma_minus_b.getFieldForLevel(level) *= masks[level];
    Ma_minus_b.toReal();


    // T^+ Mbar T n_alpha 
    // Note: we use an alias here for the sake of readability
    auto & Mbar_n_alpha = n_alpha;
    Mbar_n_alpha.toFourier();
    Mbar_n_alpha.applyPowerSpectrumFor(particle::species::all);
    Mbar_n_alpha.combineGrids();
    Mbar_n_alpha.toReal();
    for (auto level = 0; level<Nlevel; ++level)
      Mbar_n_alpha.getFieldForLevel(level) *= masksCompl[level];
    Mbar_n_alpha.toReal();

    // f = b + M (a-b) + alpha;
    auto outputs = noiseB;
    outputs.toFourier();
    outputs.applyPowerSpectrumFor(particle::species::all);
    outputs.combineGrids();
    outputs.toReal();

    outputs += Ma_minus_b;
    outputs += Mbar_n_alpha;

    // Convert back to whitenoise
    outputs.toFourier();
    outputs.applyPowerSpectrumFor(particle::species::whitenoise);
    return outputs;
  }
}

#endif //IC_SPLICE_HPP
