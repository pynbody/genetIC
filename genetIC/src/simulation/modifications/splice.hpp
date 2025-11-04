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

  // //! Return the field f which satisfies f = a in flagged region while minimising (f-b).C^-1.(f-b) elsewhere
  // template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  // fields::Field<DataType,T> spliceOneLevel(fields::Field<DataType,T> & a,
  //                                          fields::Field<DataType,T> & b,
  //                                          const fields::Field<DataType,T> & cov) {

  //     // To understand the implementation below, first read Appendix A of Cadiou et al (2021),
  //     // and/or look at the 1D toy implementation (in tools/toy_implementation/gene_splicing.ipynb) which
  //     // contains a similar derivation and near-identical implementation.

  //     assert (&a.getGrid() == &b.getGrid());
  //     assert (&a.getGrid() == &cov.getGrid());
  //     auto mask = generateMaskFromFlags(a.getGrid());
  //     auto maskCompl = generateMaskComplementFromFlags(a.getGrid());

  //     a.toFourier();
  //     b.toFourier();

  //     // The preconditioner should be almost equal to the covariance.
  //     // We however set the fundamental of the power spectrum to a non-null value,
  //     // otherwise, the mean value in the spliced region is unconstrained.
  //     fields::Field<DataType,T> preconditioner(cov);
  //     preconditioner.setFourierCoefficient(0, 0, 0, 1);

  //     fields::Field<DataType,T> delta_diff = b-a;
  //     delta_diff.applyTransferFunction(preconditioner, 0.5);
  //     delta_diff.toReal();


  //     fields::Field<DataType,T> z(delta_diff);
  //     z*=mask;
  //     z.toFourier();
  //     z.applyTransferFunction(preconditioner, -1.0);
  //     z.toReal();
  //     z*=maskCompl;
  //     z.toFourier();
  //     z.applyTransferFunction(preconditioner, 0.5);
  //     z.toReal();

  //     auto X = [&preconditioner, &maskCompl](const fields::Field<DataType,T> & input) -> fields::Field<DataType,T>
  //     {
  //       fields::Field<DataType,T> v(input);
  //       assert (!input.isFourier());
  //       assert (!v.isFourier());
  //       v.toFourier();
  //       v.applyTransferFunction(preconditioner, 0.5);
  //       v.toReal();
  //       v*=maskCompl;
  //       v.toFourier();
  //       v.applyTransferFunction(preconditioner, -1.0);
  //       v.toReal();
  //       v*=maskCompl;
  //       v.toFourier();
  //       v.applyTransferFunction(preconditioner, 0.5);
  //       v.toReal();
  //       return v;
  //     };


  //     fields::Field<DataType,T> alpha = tools::numerics::conjugateGradient<DataType>(X,z);

  //     alpha.toFourier();
  //     alpha.applyTransferFunction(preconditioner, 0.5);
  //     alpha.toReal();

  //     fields::Field<DataType,T> bInDeltaBasis(b);
  //     bInDeltaBasis.toFourier();
  //     bInDeltaBasis.applyTransferFunction(preconditioner, 0.5);
  //     bInDeltaBasis.toReal();

  //     alpha*=maskCompl;
  //     alpha+=bInDeltaBasis;

  //     delta_diff*=mask;
  //     alpha-=delta_diff;

  //     assert(!alpha.isFourier());
  //     alpha.toFourier();
  //     alpha.applyTransferFunction(preconditioner, -0.5);
  //     alpha.toReal();

  //     return alpha;
  // }

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

      fields::OutputField<T> delta_diff(b);
      delta_diff -= a;
      for (size_t level=0; level<Nlevel; ++level) {
        auto & input = delta_diff.getFieldForLevel(level);
        input.applyTransferFunction(covs[level], 0.5);
        input.toReal();
      }

      fields::OutputField<T> z(delta_diff);
      for (size_t level=0; level<Nlevel; ++level) {
        auto & z_level = z.getFieldForLevel(level);
        z_level*=masks[level];
        z_level.toFourier();
        z_level.applyTransferFunction(covs[level], -1.0);
        z_level.toReal();
        z_level*=masksCompl[level];
        z_level.toFourier();
        z_level.applyTransferFunction(covs[level], 0.5);
        z_level.toReal();
      }

      auto X = [&](const fields::OutputField<T> & inputs) -> fields::OutputField<T>
      {
        auto outputs = inputs;
        for (size_t level=0; level<Nlevel; ++level) {
          auto& v = outputs.getFieldForLevel(level);
          assert (!v.isFourier());
          v.toFourier();
          v.applyTransferFunction(covs[level], 0.5);
          v.toReal();
          v*=masksCompl[level];
          v.toFourier();
          v.applyTransferFunction(covs[level], -1.0);
          v.toReal();
          v*=masksCompl[level];
          v.toFourier();
          v.applyTransferFunction(covs[level], 0.5);
          v.toReal();
        }
        return outputs;
      };


      fields::OutputField<DataType> alpha = tools::numerics::conjugateGradient<DataType>(X,z);

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
