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

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> Mbar_Cm1_Mbar(
    const fields::OutputField<DataType> & inputs,
    const auto& covs, const auto& masks, const auto& masksCompl, size_t Nlevel
  ) {
    auto outputs = inputs;
    // outputs = T_op_T(outputs, [&](const int level, fields::Field<DataType,T> & input) {
    for (size_t level = 0; level < Nlevel; ++level) {
      auto & input = outputs.getFieldForLevel(level);
      input.toFourier();
      input.applyTransferFunction(covs[level], 0.5);
      input.toReal();
      input *= masksCompl[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= masksCompl[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], 0.5);
      input.toReal();
    }
    // });
    outputs.toReal();
    return outputs;
  }

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> Mbar_Cm1_M(
    const fields::OutputField<DataType> & inputs,
    const auto& covs, const auto& masks, const auto& masksCompl, size_t Nlevel
  ) {
    auto outputs = inputs;
    // outputs = T_op_T(outputs, [&](const int level, fields::Field<DataType,T> & input) {
    for (size_t level = 0; level < Nlevel; ++level) {
      auto & input = outputs.getFieldForLevel(level);
      input.toFourier();
      input.applyTransferFunction(covs[level], 0.5);
      input.toReal();
      input *= masks[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], -1.0);
      input.toReal();
      input *= masksCompl[level];
      input.toFourier();
      input.applyTransferFunction(covs[level], 0.5);
      input.toReal();
    // });
    }
    outputs.toReal();
    return outputs;
  };

  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> splice(fields::OutputField<DataType> & a,
                                       fields::OutputField<DataType> & b) {

      auto filters = a.getFilters();

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

      a.toFourier();
      b.toFourier();

      fields::OutputField<T> delta(b);
      delta -= a;

      fields::OutputField<T> z = Mbar_Cm1_M(delta, covs, masks, masksCompl, Nlevel);

      auto A = [&](const fields::OutputField<T> & inputs) {
        return Mbar_Cm1_Mbar(inputs, covs, masks, masksCompl, Nlevel);
      };


      fields::OutputField<DataType> alpha = tools::numerics::conjugateGradient<DataType>(A, z);

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
