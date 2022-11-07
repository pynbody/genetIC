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
    assert (noiseA.isFourierOnAllLevels());
    assert (noiseB.isFourierOnAllLevels());


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
        const TOPERATOR Top,
        const bool filter_in_window=true
    ) -> fields::OutputField<DataType> {
      fields::OutputField<DataType> outputs(inputs.getContext(), inputs.getTransferType());
      fields::OutputField<DataType> inputCopy(inputs);
      const double power_out = (Top == TDAGGER) ? 0.5 : -0.5;

      assert (inputs.isRealOnAllLevels());

      // Compute operator using each possible pair
      for (auto level=0; level<Nlevel; ++level) {
        auto & out = outputs.getFieldForLevel(level);
        out.toReal();

        for (auto source_level=std::max(0, level-1); source_level<std::min(Nlevel, level+2); ++source_level) {
          T pixel_volume_ratio = multiLevelContext.getWeightForLevel(level) /
                                 multiLevelContext.getWeightForLevel(source_level);

          fields::Field<DataType, T> tmp(out.getGrid(), false);

          // if (source_level != level) continue;
          auto & this_level_filter = filters.getFilterForLevel(level);
          auto & source_level_filter = filters.getFilterForLevel(source_level);

          if (source_level > level) {        // Contribution from finer level
            tmp.deInterpolate(inputCopy.getFieldForLevel(source_level));
            // tmp.addFieldFromDifferentGrid(inputCopy.getFieldForLevel(source_level));
            tmp.toFourier();
            tmp.applyFilter(source_level_filter);
            tmp.applyTransferFunction(covs[level], 0.5);
            tmp.toReal();
            op(level, tmp);
            tmp.toFourier();
            tmp.applyTransferFunction(covs[level], power_out);
            tmp.applyFilter(this_level_filter);
            tmp.toReal();
          } else if (source_level < level) { // Contribution from coarser level
            auto source_field = inputCopy.getFieldForLevel(source_level).copy();
            source_field->toFourier();
            source_field->applyFilter(source_level_filter);
            source_field->applyTransferFunction(covs[source_level], 0.5);
            source_field->toReal();
            op(source_level, *source_field);
            source_field->toFourier();
            source_field->applyTransferFunction(covs[source_level], power_out);
            source_field->applyFilter(this_level_filter);
            source_field->toReal();
            tmp.toReal();
            tmp.addFieldFromDifferentGrid(*source_field);
          } else {
            tmp += inputCopy.getFieldForLevel(level);
            tmp.toFourier();
            if (filter_in_window && level < Nlevel-1) {
              auto window = multiLevelContext.getGridForLevel(level+1).getWindow();
              tmp.applyFilterInWindow(this_level_filter, window, true);
            } else {
              tmp.applyFilter(this_level_filter);
            }
            tmp.applyTransferFunction(covs[level], 0.5);
            tmp.toReal();
            op(level, tmp);
            tmp.toFourier();
            tmp.applyTransferFunction(covs[level], power_out);
            if (filter_in_window && level < Nlevel-1) {
              auto window = multiLevelContext.getGridForLevel(level+1).getWindow();
              tmp.applyFilterInWindow(this_level_filter, window, false);
            } else {
              tmp.applyFilter(this_level_filter);
            }
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
      const bool precondition = false;

      if (precondition) {
        Tdagger_op_T(inputs, [covs](const int level, fields::Field<DataType, T> & input) {
          input.toFourier();
          input.applyTransferFunction(covs[level], 0.5);
        });
      }
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
    rhs.getFieldForLevel(0).dumpGridData("rhs-0.npz");
    if (Nlevel > 1) rhs.getFieldForLevel(1).dumpGridData("rhs-1.npz");
    auto n_alpha = tools::numerics::minres<DataType>(Q, rhs);
    // auto n_alpha = tools::numerics::bicgstab<DataType>(Q, rhs);
    preconditioner(n_alpha);

    // if (Nlevel > 1) {
    //   auto Ntot = rhs.getFieldForLevel(0).getGrid().size3 + ((Nlevel > 1) ? rhs.getFieldForLevel(1).getGrid().size3 : 0);
    //   for (auto i = 0; i < Ntot; ++i) {
    //     if (i % 10 == 0) std::cout << i << "/" << Ntot << std::endl;
    //     auto input = rhs;
    //     input *= 0;
    //     if (i < rhs.getFieldForLevel(0).getGrid().size3)
    //       input.getFieldForLevel(0)[i] = 1;
    //     else
    //       input.getFieldForLevel(1)[i-rhs.getFieldForLevel(0).getGrid().size3] = 1;
    //     auto output = Q(input);
    //     std::ostringstream stream;
    //     stream << "field" << i;
    //     std::string filename = stream.str();
    //     output.getFieldForLevel(0).dumpGridData(filename + "_0.npz");
    //     output.getFieldForLevel(1).dumpGridData(filename + "_1.npz");
    //   }
    // }

    // T^+ M T (n_a-n_b)
    auto Ma_minus_b = noiseA;
    Ma_minus_b -= noiseB;
    Ma_minus_b.toReal();
    Ma_minus_b = Tplus_op_T(Ma_minus_b, [&masks](const int level, fields::Field<DataType, T> & input) {
      input.toReal();
      input *= masks[level];
    });

    // T^+ Mbar T n_alpha
    auto Mbar_n_alpha = Tplus_op_T(n_alpha, [&masksCompl](const int level, fields::Field<DataType, T> & input) {
      input.toReal();
      input *= masksCompl[level];
    });

    // Reconstruct solution
    auto outputs = noiseB;
    // outputs.toReal();
    // outputs = Tplus_op_T(outputs, [](const int level, fields::Field<DataType, T> & input) {});

    Ma_minus_b.toReal();
    Mbar_n_alpha.toReal();
    outputs.toReal();
    for (auto level = 0; level < Nlevel; ++level) {
      auto & field = outputs.getFieldForLevel(level);
      field += Ma_minus_b.getFieldForLevel(level);
      field += Mbar_n_alpha.getFieldForLevel(level);
    }

    outputs.toFourier();
    return outputs;
  }
}

#endif //IC_SPLICE_HPP
