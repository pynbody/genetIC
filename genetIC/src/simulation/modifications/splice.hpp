#ifndef IC_SPLICE_HPP
#define IC_SPLICE_HPP

#include <complex>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/numerics/cg.hpp>

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


      fields::Field<DataType,T> alpha = tools::numerics::conjugateGradient<DataType>(X,z);

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
    std::vector<fields::Field<char, T>> Mbars;
    std::vector<fields::Field<char, T>> Ms;
    std::vector<fields::Field<DataType, T>> covs;

    auto transferType = noiseA.getTransferType();
    auto filters = noiseA.getFilters();
    
    for(size_t level=0; level<Nlevel; ++level) {
      Mbars.push_back(generateMaskComplementFromFlags(multiLevelContext.getGridForLevel(level)));
      Ms.push_back(generateMaskFromFlags(multiLevelContext.getGridForLevel(level)));

      fields::Field<DataType, T> cov(*multiLevelContext.getCovariance(level, particle::species::all));
      cov.setFourierCoefficient(0, 0, 0, 1);
      covs.push_back(cov);
    }

    // Compute the operator that applies T^+ ... T *at all levels*
    auto Tplus_op_T = [&covs, &filters, &multiLevelContext, &transferType, Nlevel](
        fields::OutputField<DataType> &inputs,
        const std::function<void(const int, fields::Field<DataType, T> &)> op
    ) -> fields::OutputField<DataType> {
      fields::OutputField<DataType> outputs(inputs);
      
      // TODO: optimize, as we don't need to compute the off-diagonal term twice (operation is symmetric)
      // TODO: precompute F C^0.5 for each level (otherwise it is recomputed multiple times)

      // Compute operator using each possible pair
      for (auto level=0; level<Nlevel; ++level) {
        auto & cov_me = covs[level];
        auto & f_me = filters.getFilterForLevel(level);
        auto & out = outputs.getFieldForLevel(level);

        fields::Field<DataType, T> tmp2(multiLevelContext.getGridForLevel(level), false);

        // This loop computes contribution from off-diagonal terms (high-from-low and low-from-high)
        // as well as the diagonal term (low-from-low).
        for (auto level_other=0; level_other<Nlevel; ++level_other) {
          T pixel_volume_ratio = multiLevelContext.getWeightForLevel(level) /
                                 multiLevelContext.getWeightForLevel(level_other);
          auto & f_other = filters.getFilterForLevel(level_other);
          auto & cov_other = covs[level_other];
          auto tmp_from_other = inputs.getFieldForLevel(level_other).copy();

          tmp_from_other->toFourier();
          tmp_from_other->applyFilter(f_other);
          tmp_from_other->applyTransferFunction(cov_other, 0.5);
          tmp_from_other->toReal();

          tmp2 *= 0;
          if (level_other < level) {             // Contribution from coarser levels
            tmp2.addFieldFromDifferentGrid(*tmp_from_other);  // upsample
            op(level, tmp2);                                  // then apply operator
          } else if (level_other > level) {      // Contribution from finer levels
            op(level_other, *tmp_from_other);                 // apply operator
            tmp_from_other->toReal();
            tmp2.addFieldFromDifferentGrid(*tmp_from_other);  // then downsample
          } else if (level_other == level) {     // Contribution from same level
            tmp2 += *tmp_from_other;
            op(level, tmp2);
          }

          tmp2.toFourier();
          tmp2.applyTransferFunction(cov_me, 0.5);
          tmp2.applyFilter(f_me);

          tmp2.toReal();
          out.addScaled(tmp2, sqrt(pixel_volume_ratio));
        }
      }

      return outputs;
    };

    // Apply C^0.5 *at a single level*
    auto Chalf_op = [&covs, &multiLevelContext](const int level, fields::Field<DataType,T> & input) {
      input.toFourier();
      input.applyTransferFunction(covs[level], 0.5);
    };

    // Apply Mbar C^-1 Mbar *at a single level*
    auto Mbar_Cm1_Mbar_op = [&covs, &multiLevelContext, &Mbars](const int level, fields::Field<DataType,T> & input) {
      input.toReal();
      input *= Mbars[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= Mbars[level];
    };

    // Apply Mbar C^-1 M *at a single field*
    auto Mbar_Cm1_M_op = [&covs, &multiLevelContext, &Mbars, &Ms](const int level, fields::Field<DataType,T> & input) {
      input.toReal();
      input *= Ms[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= Mbars[level];
    };

    // Full splicing operator
    // Q := T^+ Mbar C^-1 Mbar T (+ preconditioning)
    auto Q = [&Tplus_op_T, &Mbar_Cm1_Mbar_op, &Chalf_op](
      fields::OutputField<DataType> & inputs
    ) -> fields::OutputField<DataType> {
      // Precondition
      auto outputs = Tplus_op_T(inputs, Chalf_op);
      outputs      = Tplus_op_T(outputs, Mbar_Cm1_Mbar_op);
      outputs      = Tplus_op_T(outputs, Chalf_op);
      outputs.toReal();
      return outputs;
    };

    // Compute the rhs T^+ Mbar C^-1 M T (b - a)
    fields::OutputField<T> rhs(noiseB);
    rhs.addScaled(noiseA, -1);
    Tplus_op_T(rhs, Mbar_Cm1_M_op);
    Tplus_op_T(rhs, Chalf_op); // + preconditioning

    // Solve the linear system [T^+ Mbar C^-1 Mbar T] nZ = rhs
    rhs.toReal();
    auto nZ = tools::numerics::conjugateGradient2<DataType>(Q, rhs);
    Tplus_op_T(nZ, Chalf_op); // + preconditioning

    return nZ;
    // Compute the solution
  }
}

#endif //IC_SPLICE_HPP
