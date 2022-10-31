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


      logging::entry() << "Solving splicing problem ------------------------" << std::endl;
      fields::Field<DataType,T> alpha = tools::numerics::minresField<DataType>(X, z);

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
      covs.push_back(cov);
    }

    // Compute the operator that applies T^+ ... T *at all levels*
    auto T_op_T = [&covs, &filters, &multiLevelContext, Nlevel](
        const fields::OutputField<DataType> &inputs,
        const std::function<fields::Field<DataType, T>(const int, const fields::Field<DataType, T> &)> op,
        const TOPERATOR Top
    ) -> fields::OutputField<DataType> {
      fields::OutputField<DataType> outputs(inputs.getContext(), inputs.getTransferType());
      fields::OutputField<DataType> inputCopy(inputs);
      
      // TODO: optimize, as we don't need to compute the off-diagonal term twice (operation is symmetric)
      // TODO: precompute F C^0.5 for each level (otherwise it is recomputed multiple times)

      // Compute operator using each possible pair
      for (auto level=0; level<Nlevel; ++level) {
        auto & cov_me = covs[level];
        auto & f_me = filters.getFilterForLevel(level);
        auto & out = outputs.getFieldForLevel(level);
        out.toReal();

        fields::Field<DataType, T> tmp(multiLevelContext.getGridForLevel(level), false);

        // This loop computes contribution from off-diagonal terms (high-from-low and low-from-high)
        // as well as the diagonal term (low-from-low).
        for (auto level_other=0; level_other<Nlevel; ++level_other) {
          T pixel_volume_ratio = multiLevelContext.getWeightForLevel(level) /
                                 multiLevelContext.getWeightForLevel(level_other);
          auto & f_other = filters.getFilterForLevel(level_other);
          auto & cov_other = covs[level_other];
          auto tmp_from_other = inputCopy.getFieldForLevel(level_other).copy();

          tmp_from_other->toFourier();
          std::cout << " " << f_other << " "<< " C^0.5";
          tmp_from_other->applyFilter(f_other);
          tmp_from_other->applyTransferFunction(cov_other, 0.5);
          tmp_from_other->toReal();

          tmp *= 0;
          if (level_other < level) {             // Contribution from coarser level
            // Upsample then compute operation
            std::cout << " P_" << level_other;
            tmp.addFieldFromDifferentGrid(*tmp_from_other); 
            tmp = op(level, tmp);
          } else if (level_other > level) {      // Contribution from finer level
            // Compute operation then downsample
            tmp.addFieldFromDifferentGrid(op(level_other, *tmp_from_other));
            std::cout << " P^+_" << level_other;
          } else if (level_other == level) {     // Contribution from same level
            tmp = op(level, *tmp_from_other);
          }

          tmp.toFourier();
          if (Top == TDAGGER) {
            std::cout << " C^0.5 ";
            tmp.applyTransferFunction(cov_me, 0.5);
          } else if (Top == TPLUS) {
            std::cout << " C^-0.5 ";
            tmp.applyTransferFunction(cov_me, -0.5);
          }
          std::cout << f_me;
          tmp.applyFilter(f_me);

          tmp.toReal();
          std::cout << " * " << sqrt(pixel_volume_ratio);
          out.addScaled(tmp, sqrt(pixel_volume_ratio));
        }
      }

      return outputs;
    };

    auto Tplus_op_T = [&T_op_T](const fields::OutputField<DataType> &inputs,
                                const std::function<fields::Field<DataType, T>(const int, const fields::Field<DataType, T> &)> op
    ) -> fields::OutputField<DataType> {
      return T_op_T(inputs, op, TPLUS);
    };

    auto Tdagger_op_T = [&T_op_T](const fields::OutputField<DataType> &inputs,
                                  const std::function<fields::Field<DataType, T>(const int, const fields::Field<DataType, T> &)> op
    ) -> fields::OutputField<DataType> {
      return T_op_T(inputs, op, TDAGGER);
    };

    // Compute Mbar C^-1 Mbar for a single level
    auto Mbar_Cm1_Mbar_op = [&covs, &multiLevelContext, &masksCompl](const int level, const fields::Field<DataType,T> & input) -> fields::Field<DataType,T> {
      auto output = input;
      output.toReal();
      std::cout << " Mbar";
      output *= masksCompl[level];
      output.toFourier();
      std::cout << " C^-1";
      output.applyTransferFunction(covs[level], -1.0);
      output.toReal();
      std::cout << " Mbar";
      output *= masksCompl[level];
      return output;
    };

    // Compute Mbar C^-1 M for a single field
    auto Mbar_Cm1_M_op = [&covs, &multiLevelContext, &masksCompl, &masks](const int level, const fields::Field<DataType,T> & input) -> fields::Field<DataType,T> {
      auto output = input;
      output.toReal();
      std::cout << " M";
      output *= masks[level];
      output.toFourier();
      std::cout << " C^-1";
      output.applyTransferFunction(covs[level], -1.0);
      output.toReal();
      std::cout << " Mbar";
      output *= masksCompl[level];
      return output;
    };

    // Full splicing operator
    // Q := T^+ Mbar C^-1 Mbar T
    auto Q = [&Tdagger_op_T, &Mbar_Cm1_Mbar_op](
      const fields::OutputField<DataType> & inputs
    ) -> fields::OutputField<DataType> {
      // Precondition
      std::cout << "                                               Q = ";
      auto outputs = inputs;
      outputs = Tdagger_op_T(outputs, Mbar_Cm1_Mbar_op);
      outputs.toReal();
      std::cout << std::endl;
      return outputs;
    };

    // Compute the rhs T^+ Mbar C^-1 M T (b - a)
    fields::OutputField<T> rhs(noiseB);
    rhs.addScaled(noiseA, -1);
    rhs.toReal();


    std::cout << "                                               nB-nA";
    rhs = Tdagger_op_T(rhs, Mbar_Cm1_M_op);
    std::cout << std::endl;
    // Solve the linear system [T^t Mbar C^-1 Mbar T] n_alpha = rhs
    // that approximates [T^+ Mbar C^-1 Mbar T] and is symmetric
    rhs.toReal();

    auto n_alpha = tools::numerics::minres<DataType>(Q, rhs);

    // T^+ M T (n_a-n_b)
    fields::OutputField<DataType> Ma_minus_b(noiseA);
    Ma_minus_b.addScaled(noiseB, -1);
    Ma_minus_b.applyPowerSpectrumFor(particle::species::all);
    Ma_minus_b.toReal();
    for (auto level = 0; level<Nlevel; ++level)
      Ma_minus_b.getFieldForLevel(level) *= masks[level];

    // T^+ Mbar T n_alpha
    // Note: we use an alias here for the sake of readability
    fields::OutputField<DataType> & Mbar_n_alpha = n_alpha;
    Mbar_n_alpha.applyPowerSpectrumFor(particle::species::all);
    Mbar_n_alpha.toReal();
    for (auto level = 0; level<Nlevel; ++level)
      Mbar_n_alpha.getFieldForLevel(level) *= masksCompl[level];

    // f = b + M (a-b) + alpha; 
    fields::OutputField<DataType> outputs(noiseB);
    outputs.applyPowerSpectrumFor(particle::species::all);
    outputs.toReal();

    outputs += Ma_minus_b;
    outputs += Mbar_n_alpha;

    // outputs.combineGrids();

    Ma_minus_b.getFieldForLevel(0).dumpGridData("Ma_minus_b-0.npz");
    Mbar_n_alpha.getFieldForLevel(0).dumpGridData("Mbar_n_alpha-0.npz");
    outputs.getFieldForLevel(0).dumpGridData("output-0.npz");

    // Ma_minus_b.getFieldForLevel(1).dumpGridData("Ma_minus_b-1.npz");
    // Mbar_n_alpha.getFieldForLevel(1).dumpGridData("Mbar_n_alpha-1.npz");
    // outputs.getFieldForLevel(1).dumpGridData("output-1.npz");

    // Convert back to whitenoise
    outputs.applyPowerSpectrumFor(particle::species::whitenoise);

    return outputs;
  }
}

#endif //IC_SPLICE_HPP
