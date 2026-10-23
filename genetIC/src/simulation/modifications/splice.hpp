#ifndef IC_SPLICE_HPP
#define IC_SPLICE_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <src/simulation/field/multilevelfield.hpp>
#include <src/tools/numerics/cg.hpp>
#include <stdexcept>
#include <type_traits>

namespace modifications {

  /*! Production combine and its exact transpose, following manual_splicing.ipynb
   * choices 0/1. UpSample is production cubic interpolation J; its transpose
   * scatters the same weights (not the notebook's block average DownSample).
   * OutputField stores levels from coarse to fine; latents have unit real variance.
   * the covariance fields already include the physical cell-volume normalization.
   * No context flags, power-spectrum zero modes, or input arrays are changed.
   */
  template <typename DataType> class SpliceSynthesis {
    using T = tools::datatypes::strip_complex<DataType>;
    using State = fields::OutputField<DataType>;
    const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext;
    filters::FilterFamily<T> filters;
    const std::vector<fields::Field<DataType, T>> *covs;

    std::vector<std::vector<T>> Chalf, FH, FL;

    // ======================================================================
    // Fourier multipliers: C^(1/2), high-pass FH, and low-pass FL
    // ======================================================================
    template <typename Fn> std::vector<T> makeKernel(size_t level, Fn fn) const {
      using Complex = tools::datatypes::ensure_complex<DataType>;
      fields::Field<DataType, T> coefficients(multiLevelContext.getGridForLevel(level), true);
      coefficients.forEachFourierCellInt(
          [&](Complex, int x, int y, int z) -> Complex { return Complex(fn(x, y, z), 0); });
      coefficients.ensureFourierModesAreMirrored();
      const auto &data = coefficients.getDataVector();
      std::vector<T> kernel(data.size() / 2);
#pragma omp parallel for schedule(static)
      for (size_t j = 0; j < kernel.size(); ++j)
        kernel[j] = data[2 * j];
      return kernel;
    }

    std::vector<T> makeFilterKernel(size_t level, const filters::Filter<T> &filter) const {
      const T kMin = multiLevelContext.getGridForLevel(level).getFourierKmin();
      return makeKernel(level, [&](int x, int y, int z) {
        const T kx = x * kMin, ky = y * kMin, kz = z * kMin;
        const T k = std::sqrt(double(kx * kx + ky * ky + kz * kz));
        return filter(k);
      });
    }

    // Immutable real multipliers on FFTW's interleaved Fourier storage. Apply
    // two factors in their original order, without reassociating the product.
    void applyKernels(fields::Field<DataType, T> &field, const std::vector<T> &first,
                      const std::vector<T> *second = nullptr) const {
      field.toFourier();
      auto &data = field.getDataVector();
      assert(data.size() == 2 * first.size());
      if (second) {
        assert(second->size() == first.size());
#pragma omp parallel for schedule(static)
        for (size_t j = 0; j < first.size(); ++j) {
          data[2 * j] *= first[j];
          data[2 * j + 1] *= first[j];
          data[2 * j] *= (*second)[j];
          data[2 * j + 1] *= (*second)[j];
        }
      } else {
#pragma omp parallel for schedule(static)
        for (size_t j = 0; j < first.size(); ++j) {
          data[2 * j] *= first[j];
          data[2 * j + 1] *= first[j];
        }
      }
    }

    void applyChalfFH(fields::Field<DataType, T> &field, size_t level, bool highPassFirst = false) const {
      field.toFourier();
      if (!Chalf.empty()) {
        if (!level)
          applyKernels(field, Chalf[level]);
        else if (highPassFirst)
          applyKernels(field, FH[level], &Chalf[level]);
        else
          applyKernels(field, Chalf[level], &FH[level]);
        return;
      }
      if (level && highPassFirst)
        field.applyFilter(filters.getHighPassFilterForLevel(level));
      field.applyTransferFunction(
          covs ? covs->at(level) : *multiLevelContext.getCovariance(level, particle::species::all), .5);
      if (level && !highPassFirst)
        field.applyFilter(filters.getHighPassFilterForLevel(level));
    }

    void applyFL(fields::Field<DataType, T> &field, size_t parentLevel) const {
      if (!FL.empty())
        applyKernels(field, FL[parentLevel]);
      else
        field.applyFilter(filters.getLowPassFilterForLevel(parentLevel));
    }

  public:
    // manual_splicing.ipynb: choice 1 filters before UpSample; choice 0 after.
    // Match the particle generator's build-time convention; no window choices.
#ifdef FILTER_ON_COARSE_GRID
    static constexpr bool filter_on_coarse_grid = true;
#else
    static constexpr bool filter_on_coarse_grid = false;
#endif

    explicit SpliceSynthesis(const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext,
                             const std::vector<fields::Field<DataType, T>> *covs = nullptr)
        : multiLevelContext(multiLevelContext), filters(multiLevelContext), covs(covs) {
      if (covs && covs->size() != multiLevelContext.getNumLevels())
        throw std::runtime_error("combine requires one covariance per grid level");
      if (covs)
        for (size_t level = 0; level < covs->size(); ++level)
          if (&covs->at(level).getGrid() != &multiLevelContext.getGridForLevel(level))
            throw std::runtime_error("combine covariance grids do not match the input hierarchy");
      if (multiLevelContext.getLevelsAreCombined())
        throw std::runtime_error("Splicing requires uncombined white-noise fields");
      for (size_t level = 1; level < multiLevelContext.getNumLevels(); ++level)
        if (multiLevelContext.getGridForLevel(level).cellSize >
            multiLevelContext.getGridForLevel(level - 1).cellSize)
          throw std::runtime_error("Splicing requires nondecreasing grid resolution");
    }

    // Build once per solve, after validating a nonempty selection. Other users
    // retain the original production kernels as an independent reference path.
    void prepareFourierKernels() {
      if constexpr (std::is_same<DataType, T>::value) {
        const size_t levels = multiLevelContext.getNumLevels();
        Chalf.resize(levels);
        FH.resize(levels);
        FL.resize(levels ? levels - 1 : 0);
        for (size_t level = 0; level < levels; ++level) {
          const auto generatedCov = covs ? std::shared_ptr<fields::Field<DataType, T>>()
                                         : multiLevelContext.getCovariance(level, particle::species::all);
          const auto &C = covs ? covs->at(level) : *generatedCov;
          Chalf[level] = makeKernel(level, [&](int x, int y, int z) {
            auto value = C.getFourierCoefficient(x, y, z).real();
            if (value != 0)
              value = std::sqrt(value);
            return value;
          });
          if (level)
            FH[level] = makeFilterKernel(level, filters.getHighPassFilterForLevel(level));
          if (level + 1 < levels)
            FL[level] = makeFilterKernel(filter_on_coarse_grid ? level : level + 1,
                                         filters.getLowPassFilterForLevel(level));
        }
      }
    }

    // ======================================================================
    // Forward synthesis: white noise -> combined density, coarse to fine
    // ======================================================================
    // With a selection hint, only those finest-level samples are guaranteed;
    // ancestors stay complete because descendants depend on their full fields.
    State combine(const State &a, const std::vector<size_t> *finestSelection = nullptr) const {
      State outputs(a);
      for (size_t level = 0; level < outputs.getNumLevels(); ++level) {
        auto &field = outputs.getFieldForLevel(level);
        applyChalfFH(field, level);
        if (level) {
          if constexpr (filter_on_coarse_grid) {
            // choice 1: UpSample(FL_coarse(coarse)).
            auto parent = outputs.getFieldForLevel(level - 1).copy();
            applyFL(*parent, level - 1);
            parent->toReal();
            field.addFieldFromDifferentGrid(*parent,
                                            level + 1 == outputs.getNumLevels() ? finestSelection : nullptr);
          } else {
            // choice 0: FL_fine(UpSample(coarse)).
            auto &parent = outputs.getFieldForLevel(level - 1);
            parent.toReal();
            auto inherited =
                std::make_shared<fields::Field<DataType, T>>(multiLevelContext.getGridForLevel(level), false);
            inherited->addFieldFromDifferentGrid(parent);
            applyFL(*inherited, level - 1);
            field.toFourier();
            field += *inherited;
          }
        }
      }
      outputs.toReal();
      return outputs;
    }

    // ======================================================================
    // Transpose synthesis: density -> white noise, fine to coarse
    // ======================================================================
    // Gram applications can retain Fourier white noise, avoiding an inverse
    // transform here followed immediately by a forward FFT in combine().
    State combineTranspose(const State &v, bool realOutput = true) const {
      State inputs(v);
      inputs.toReal();
      State outputs(multiLevelContext, particle::species::whitenoise);
      outputs.getFieldForLevel(0);
      for (size_t level = v.getNumLevels(); level-- > 0;) {
        auto &source = inputs.getFieldForLevel(level);
        if (level) {
          // J is interpolation of the parent at each child centroid. J^T
          // scatters with the identical cubic (or linear) weights and boundaries.
          auto scatter = std::make_shared<fields::Field<DataType, T>>(
              multiLevelContext.getGridForLevel(level - 1), false);
          // With fine-grid filtering, (FL J)^T = J^T FL. Otherwise
          // (J FL)^T = FL J^T, so filter the scatter below instead.
          std::shared_ptr<fields::Field<DataType, T>> inherited;
          const auto *toScatter = &source;
          if constexpr (!filter_on_coarse_grid) {
            inherited = source.copy();
            applyFL(*inherited, level - 1);
            inherited->toReal();
            toScatter = inherited.get();
          }
          const auto &grid = source.getGrid();
          for (size_t i = 0; i < grid.size3; ++i) {
            const auto value = (*toScatter)[i];
            // Masked dual fields are sparse. Preserve the order of nonzero
            // additions; parallelizing this overlapping scatter would race.
            if (value != DataType(0))
              scatter->deInterpolate(grid.getCentroidFromIndex(i), value);
          }
          if constexpr (filter_on_coarse_grid)
            applyFL(*scatter, level - 1);
          scatter->toReal();
          inputs.getFieldForLevel(level - 1) += *scatter;
        }
        auto &out = outputs.getFieldForLevel(level);
        out = fields::Field<DataType, T>(source);
        applyChalfFH(out, level, true);
        if (realOutput)
          out.toReal();
      }
      return outputs;
    }
  };

  template <typename DataType, typename T = tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType> combine(const fields::OutputField<DataType> &a,
                                        const std::vector<fields::Field<DataType, T>> &covs) {
    SpliceSynthesis<DataType> synthesis(a.getContext(), &covs);
    return synthesis.combine(a);
  }

  // Select flagged cells that have particles in the output. The mapper's zoom
  // masks describe particle ownership; child FFT windows are larger than this.
  // Grid-only callers instead select the finest available cell at each point.
  template <typename DataType>
  std::vector<std::vector<size_t>>
  getSpliceMask(const multilevelgrid::MultiLevelGrid<DataType> &multiLevelContext,
                const std::vector<std::vector<size_t>> *zoomMasks, int borderSafety) {
    const size_t levels = multiLevelContext.getNumLevels();
    std::vector<std::vector<size_t>> replaced;
    if (zoomMasks) {
      if (zoomMasks->size() + 1 != levels)
        throw std::runtime_error("Splice zoom ownership masks do not match the hierarchy");
      replaced = *zoomMasks;
      for (auto &cells : replaced)
        std::sort(cells.begin(), cells.end());
    }
    auto isReplaced = [&](size_t level, const auto &position) {
      const auto i = multiLevelContext.getOutputGridForLevel(level).getIndexFromPoint(position);
      return std::binary_search(replaced[level].begin(), replaced[level].end(), i);
    };
    auto hasParticle = [&](size_t level, const auto &position) {
      if (!zoomMasks)
        return level + 1 == levels || !multiLevelContext.getGridForLevel(level + 1).containsPoint(position);
      return (level == 0 || isReplaced(level - 1, position)) &&
             (level + 1 == levels || !isReplaced(level, position));
    };

    std::vector<std::vector<size_t>> selected(levels);
    for (size_t level = 0; level < levels; ++level) {
      const auto &grid = multiLevelContext.getGridForLevel(level);
      std::vector<size_t> flags;
      grid.getFlaggedCells(flags);
      for (auto i : flags) {
        const auto position = grid.getCentroidFromIndex(i);
        if (!hasParticle(level, position))
          continue;
        for (size_t ancestor = 1; ancestor <= level; ++ancestor)
          if (!multiLevelContext.getGridForLevel(ancestor).containsPointWithBorderSafety(position,
                                                                                         borderSafety))
            throw std::runtime_error("Splice region is too close to a zoom edge; enlarge the zoom grid "
                                     "(insufficient configured padding on an enclosing level)");
        selected[level].push_back(i);
      }
    }
    return selected;
  }

  /*! Minimize ||nf - nb||^2 subject to MT(nf) = MT(na), as in the notebook.
   * T is density synthesis (combine), MT = M T, and TdagMdag = T^dagger M^dagger.
   * G is the Gram operator G = (M T) (M T)^dagger, hence positive semidefinite.
   * Return nf in white-noise coordinates: the particle-output pipeline performs
   * combine afterwards. Returning the notebook's combined density here would
   * apply the transfer functions twice.
   */
  template <typename DataType, typename T = tools::datatypes::strip_complex<DataType>>
  fields::OutputField<DataType>
  splice(const fields::OutputField<DataType> &na, const fields::OutputField<DataType> &nb, T rtol,
         const std::vector<std::vector<size_t>> *zoomMasks = nullptr, int borderSafety = 4) {
    using State = fields::OutputField<DataType>;

    // ======================================================================
    // 1. Check inputs and initialise the density synthesis
    // ======================================================================
    if (!std::isfinite(rtol) || rtol <= 0 || rtol >= 1)
      throw std::runtime_error("Splice relative tolerance must be finite and between zero and one");
    if (!na.isCompatible(nb) || na.getTransferType() != particle::species::whitenoise ||
        nb.getTransferType() != particle::species::whitenoise)
      throw std::runtime_error("Splicing requires compatible white-noise fields");

    if (borderSafety < 1)
      throw std::runtime_error("Splice border safety must be positive");
    // The command-line path retains four cells. A smaller explicit value is
    // available to numerical stress tests of nearly window-filling patches.
    const auto &multiLevelContext = na.getContext();
    SpliceSynthesis<DataType> synthesis(multiLevelContext);
    // ======================================================================
    // 2. Build the diagonal selection matrix M (stored as indices and bytes)
    // ======================================================================
    const auto selected = getSpliceMask(multiLevelContext, zoomMasks, borderSafety);
    size_t n = 0;
    for (const auto &cells : selected)
      n += cells.size();
    if (!n)
      return State(nb);
    synthesis.prepareFourierKernels();

    // Cache membership once. Mask in place instead of copying every full field
    // on every operator application; keep selected values and clear FFT padding.
    std::vector<std::vector<unsigned char>> membership(selected.size());
    for (size_t level = 0; level < selected.size(); ++level) {
      membership[level].resize(multiLevelContext.getGridForLevel(level).size3, 0);
      for (auto i : selected[level])
        membership[level][i] = 1;
    }
    // ======================================================================
    // 3. Define MT = M T and TdagMdag = T^dagger M^dagger
    // ======================================================================
    auto M = [&](State &state) {
      state.toReal();
      for (size_t level = 0; level < state.getNumLevels(); ++level) {
        auto &field = state.getFieldForLevel(level);
        auto &data = field.getDataVector();
        const auto &keep = membership[level];
#pragma omp parallel for
        for (size_t i = 0; i < data.size(); ++i)
          if (i >= keep.size() || !keep[i])
            data[i] = DataType(0);
      }
    };
    auto MT = [&](const State &state) {
      auto result = synthesis.combine(state, &selected.back());
      M(result);
      return result;
    };
    // M^dagger is implicit: v is already zero-embedded at the selected cells.
    // CG preserves that subspace, so only T^dagger remains to be applied.
    auto TdagMdag = [&](const State &v, bool realOutput = true) {
      return synthesis.combineTranspose(v, realOutput);
    };
    // G, the Gram operator G = (M T) (M T)^dagger.
    auto G = [&](const State &v) { return MT(TdagMdag(v, false)); };

    // ======================================================================
    // 4. Solve G lam = d, then return nf = nb + TdagMdag(lam) in white-noise units
    // ======================================================================
    // Evaluate d = MT(na) - MT(nb) in density units, as for the final residual.
    const State MTa = MT(na);
    State d(MTa);
    const auto MTb = MT(nb);
    d -= MTb;
    const T norm_d = tools::numerics::norm(d);
    // Absolute tolerance for FFT rounding: a repeated splice may already match
    // to machine precision, making rtol * norm_d unattainably small.
    const T atol = 8 * std::numeric_limits<T>::epsilon() *
                   std::max(tools::numerics::norm(MTa), tools::numerics::norm(MTb));
    if (norm_d <= atol)
      return State(nb);
    State nf(nb);
    auto densityResidual = [&](const State &lam) {
      nf = State(nb);
      nf += TdagMdag(lam);
      State residual(MTa);
      residual -= MT(nf);
      return residual;
    };
    // Shared CG checks this actual output-density residual before accepting.
    tools::numerics::conjugateGradient<DataType>(G, d, rtol, atol, densityResidual,
                                                 std::max(size_t(100), 4 * n));
    nf.toFourier();
    return nf;
  }
} // namespace modifications
#endif
