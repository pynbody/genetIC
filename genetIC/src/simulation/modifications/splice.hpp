#ifndef IC_SPLICE_HPP
#define IC_SPLICE_HPP

#include <src/simulation/field/multilevelfield.hpp>
#include <src/tools/numerics/cg.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace modifications {

  /*! The actual density synthesis used by MultiLevelParticleGenerator, and its
   * Euclidean transpose. Latents are unit-variance real white noise on each grid;
   * the covariance fields already include the physical cell-volume normalization.
   * No context flags, power-spectrum zero modes, or input arrays are changed.
   */
  template<typename DataType>
  class SpliceSynthesis {
    using T = tools::datatypes::strip_complex<DataType>;
    using State = fields::OutputField<DataType>;
    const multilevelgrid::MultiLevelGrid<DataType> &context;
    filters::FilterFamily<T> filters;

    void colour(fields::Field<DataType, T> &field, size_t level) const {
      field.toFourier();
      field.applyTransferFunction(*context.getCovariance(level, particle::species::all), 0.5);
    }

  public:
    explicit SpliceSynthesis(const multilevelgrid::MultiLevelGrid<DataType> &context)
      : context(context), filters(context) {
      if(context.getLevelsAreCombined())
        throw std::runtime_error("Splicing requires uncombined white-noise fields");
      for(size_t l=1; l<context.getNumLevels(); ++l)
        if(context.getGridForLevel(l).cellSize > context.getGridForLevel(l-1).cellSize)
          throw std::runtime_error("Splicing requires nondecreasing grid resolution");
    }

    State forward(const State &noise) const {
      State result(noise);
      for(size_t l=0; l<result.getNumLevels(); ++l) {
        auto &field = result.getFieldForLevel(l);
        colour(field, l);
        if(l) {
          field.applyFilter(filters.getHighPassFilterForLevel(l));
          field.addFieldFromDifferentGridWithFilter(result.getFieldForLevel(l-1),
                                                    filters.getLowPassFilterForLevel(l-1));
        }
      }
      result.toReal();
      return result;
    }

    State adjoint(const State &density) const {
      State accumulated(density);
      accumulated.toReal();
      State result(context, particle::species::whitenoise);
      result.getFieldForLevel(0);
      for(size_t l=density.getNumLevels(); l-- > 0;) {
        auto &source = accumulated.getFieldForLevel(l);
        if(l) {
          // J is interpolation of the parent at each child centroid. J^T
          // scatters with the identical cubic (or linear) weights and boundaries.
          auto scatter = std::make_shared<fields::Field<DataType,T>>(context.getGridForLevel(l-1), false);
#ifndef FILTER_ON_COARSE_GRID
          auto inherited = source.copy();
          inherited->applyFilter(filters.getLowPassFilterForLevel(l-1));
          inherited->toReal();
#endif
          const auto &grid = source.getGrid();
          for(size_t i=0; i<grid.size3; ++i) {
#ifdef FILTER_ON_COARSE_GRID
            scatter->deInterpolate(grid.getCentroidFromIndex(i), source[i]);
#else
            scatter->deInterpolate(grid.getCentroidFromIndex(i), (*inherited)[i]);
#endif
          }
#ifdef FILTER_ON_COARSE_GRID
          scatter->applyFilter(filters.getLowPassFilterForLevel(l-1));
#endif
          scatter->toReal();
          accumulated.getFieldForLevel(l-1) += *scatter;
        }
        auto &out = result.getFieldForLevel(l);
        out = fields::Field<DataType,T>(source);
        if(l) out.applyFilter(filters.getHighPassFilterForLevel(l));
        colour(out, l);
        out.toReal();
      }
      return result;
    }
  };

  /*! Minimize ||n - n_environment||^2 subject to M T n = M T n_donor.
   * D=M T (M is the selected-cell mask), so D D^T is PSD by construction, including
   * singular spectra. Keep the answer in latent coordinates: normal output
   * synthesis must still run, and must run only once.
   */
  template<typename DataType, typename T=tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> splice(const fields::OutputField<DataType> &donor,
                                      const fields::OutputField<DataType> &environment,
                                      T accuracy,
                                      const std::vector<std::vector<size_t>> *zoomMasks = nullptr,
                                      int borderSafety = 4) {
    using State = fields::OutputField<DataType>;
    if(!std::isfinite(accuracy) || accuracy <= 0 || accuracy >= 1)
      throw std::runtime_error("Splice accuracy must be finite and between zero and one");
    if(!donor.isCompatible(environment) ||
       donor.getTransferType() != particle::species::whitenoise ||
       environment.getTransferType() != particle::species::whitenoise)
      throw std::runtime_error("Splicing requires compatible white-noise fields");

    if(borderSafety < 1) throw std::runtime_error("Splice border safety must be positive");
    // The command-line path retains four cells. A smaller explicit value is
    // available to numerical stress tests of nearly window-filling patches.
    const auto &context = donor.getContext();
    SpliceSynthesis<DataType> synthesis(context);
    // The front end supplies the actual parent cells replaced by particles on
    // each child. A full child FFT window is larger than its particle region.
    std::vector<std::vector<size_t>> replaced;
    if(zoomMasks) {
      replaced = *zoomMasks;
      if(replaced.size()+1 != donor.getNumLevels())
        throw std::runtime_error("Splice zoom ownership masks do not match the hierarchy");
      for(auto &cells: replaced) std::sort(cells.begin(), cells.end());
    }
    std::vector<std::vector<size_t>> selected(donor.getNumLevels());
    size_t count = 0;
    for(size_t l=0; l<donor.getNumLevels(); ++l) {
      const auto &grid = context.getGridForLevel(l);
      std::vector<size_t> flags;
      grid.getFlaggedCells(flags);
      for(auto i: flags) {
        auto position = grid.getCentroidFromIndex(i);
        if(zoomMasks) {
          if(l) {
            const auto &parent = context.getOutputGridForLevel(l-1);
            auto parentIndex = parent.getIndexFromPoint(position);
            if(!std::binary_search(replaced[l-1].begin(), replaced[l-1].end(), parentIndex))
              continue; // This child sample has no particle in the output.
          }
          if(l+1 < donor.getNumLevels()) {
            auto outputIndex = context.getOutputGridForLevel(l).getIndexFromPoint(position);
            if(std::binary_search(replaced[l].begin(), replaced[l].end(), outputIndex))
              continue; // The coarse particle was replaced by child particles.
          }
        } else if(l+1 < donor.getNumLevels() && context.getGridForLevel(l+1).containsPoint(position)) {
          // Standalone grid-only callers have no particle mapper: use leaf cells.
          continue;
        }
        for(size_t ancestor=1; ancestor<=l; ++ancestor)
          if(!context.getGridForLevel(ancestor).containsPointWithBorderSafety(position, borderSafety))
            throw std::runtime_error("Splice region is too close to a zoom edge; enlarge the zoom grid (insufficient configured padding on an enclosing level)");
        selected[l].push_back(i);
      }
      count += selected[l].size();
    }
    if(!count) return State(environment);

    auto mask = [&](State &state) {
      state.toReal();
      for(size_t l=0; l<state.getNumLevels(); ++l) {
        auto &field = state.getFieldForLevel(l);
        auto saved = field.getDataVector();
        std::fill(field.getDataVector().begin(), field.getDataVector().end(), DataType(0));
        for(auto i: selected[l]) field[i] = saved[i];
      }
    };
    auto D = [&](const State &state) {
      auto result = synthesis.forward(state);
      mask(result);
      return result;
    };
    auto G = [&](const State &state) { return D(synthesis.adjoint(state)); };

    const State target = D(donor);
    State rhs(target);
    const auto exterior = D(environment);
    rhs -= exterior;
    const T scale = tools::numerics::norm(rhs);
    // A repeated splice can have a discrepancy consisting only of FFT rounding.
    // A relative tolerance on that discrepancy alone would ask for sub-roundoff
    // accuracy. Measure the floor in density units, not latent/noise units.
    const T roundoff = 8*std::numeric_limits<T>::epsilon() *
      std::max(tools::numerics::norm(target), tools::numerics::norm(exterior));
    if(scale <= roundoff) return State(environment);
    State answer(environment);
    auto densityResidual = [&](const State &lambda) {
      answer = State(environment);
      answer += synthesis.adjoint(lambda);
      State residual(target);
      residual -= D(answer);
      return residual;
    };
    // Shared CG checks this actual output-density residual before accepting.
    tools::numerics::conjugateGradient<DataType>(G, rhs, accuracy, roundoff,
                                                densityResidual, std::max(size_t(100), 4*count));
    answer.toFourier();
    return answer;
  }
}
#endif
