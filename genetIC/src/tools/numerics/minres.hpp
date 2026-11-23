#ifndef IC_MINRES_HPP
#define IC_MINRES_HPP

#include <complex>
#include <src/tools/logging.hpp>
#include <src/tools/numerics/fourier.hpp>
#include "cg.hpp"

namespace tools {
  namespace numerics {

    /*! Solve Qx=b by MINRES, starting from zero. Q must be symmetric, but
     * may be indefinite or singular if the system is consistent. Uses Lanczos
     * tridiagonalization and Givens rotations (Paige and Saunders, 1975):
     * https://web.stanford.edu/group/SOL/software/minres/
     *
     * The signature and stopping criterion match conjugateGradient: accept only
     * when ||b-Qx|| <= max(rtol*||b||, atol). recomputeResidual optionally
     * evaluates the same residual through the actual returned observable.
     * Zero maxIterations uses dimension+1. Inconsistent systems, breakdowns,
     * and iteration exhaustion throw rather than returning an unconverged x.
     * precondition applies an SPD approximate inverse, as in CG. Lanczos uses
     * its induced inner product, preserving symmetry; acceptance still uses
     * the Euclidean norm of the original residual.
     */
    template<typename T>
    fields::OutputField<T> minres(
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> Q,
        const fields::OutputField<T> &b,
        double rtol = 1e-6, double atol = 1e-12,
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> recomputeResidual = {},
        size_t maxIterations = 0,
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> precondition = {}) {
      using State = fields::OutputField<T>;
      if(!std::isfinite(rtol) || !std::isfinite(atol) || rtol < 0 || atol < 0)
        throw std::runtime_error("MINRES tolerances must be finite and nonnegative");

      const auto checkedNorm = [](const State &state) {
        const double squared = innerProduct(state, state);
        if(!std::isfinite(squared) || squared < 0)
          throw std::runtime_error("MINRES residual became nonfinite");
        return std::sqrt(squared);
      };
      State x(b.getContext(), b.getTransferType());
      x.getFieldForLevel(0); // Trigger allocation, as in CG.
      x.toReal(); // Inner products and Lanczos vectors use real-space fields.
      const double scale = checkedNorm(b);
      const double tolerance = std::max(rtol*scale, atol);
      if(scale <= tolerance) {
        x.toReal();
        logging::entry() << "MINRES ended after 0 iterations" << std::endl;
        return x;
      }
      if(!maxIterations) {
        size_t dimension = 0;
        for(size_t level=0; level<b.getNumLevels(); ++level)
          dimension += b.getContext().getGridForLevel(level).size3;
        maxIterations = dimension+1;
      }

      const auto applyPreconditioner = [&](const State &r) {
        return precondition ? precondition(r) : State(r);
      };
      const auto metricNorm = [&](const State &r, const State &weighted) {
        const double squared = innerProduct(r, weighted);
        if(!std::isfinite(squared) || squared < 0 || (squared == 0 && checkedNorm(r) != 0))
          throw std::runtime_error("MINRES preconditioner must be positive definite");
        return std::sqrt(squared);
      };
      State z(b), v = applyPreconditioner(b);
      const double metricScale = metricNorm(z, v);
      const double metricTolerance = (tolerance/scale)*metricScale;
      z /= metricScale;
      v /= metricScale;
      State previousZ(x), direction(x), previousDirection(x);
      double beta = 0, dbar = 0, epsilon = 0;
      double cosine = -1, sine = 0, residualEstimate = metricScale;
      size_t restarts = 0;
      for(size_t iteration=0; iteration<maxIterations; ++iteration) {
        // Three-term Lanczos recurrence for the symmetric operator.
        State nextZ = Q(v);
        nextZ.addScaled(previousZ, -beta);
        const double alpha = innerProduct(v, nextZ);
        if(!std::isfinite(alpha))
          throw std::runtime_error("MINRES Lanczos coefficient became nonfinite");
        nextZ.addScaled(z, -alpha);
        State nextV = applyPreconditioner(nextZ);
        const double nextBeta = metricNorm(nextZ, nextV);

        // Apply the preceding rotations to the new tridiagonal column.
        const double oldEpsilon = epsilon;
        const double delta = cosine*dbar + sine*alpha;
        const double gbar = sine*dbar - cosine*alpha;
        epsilon = sine*nextBeta;
        dbar = -cosine*nextBeta;
        const double gamma = std::hypot(gbar, nextBeta);
        if(!std::isfinite(gamma) || gamma == 0)
          throw std::runtime_error("MINRES breakdown before satisfying tolerance");
        cosine = gbar/gamma;
        sine = nextBeta/gamma;
        const double step = cosine*residualEstimate;
        residualEstimate *= sine;

        State nextDirection(v);
        nextDirection.addScaled(previousDirection, -oldEpsilon);
        nextDirection.addScaled(direction, -delta);
        nextDirection /= gamma;
        x.addScaled(nextDirection, step);
        previousDirection = direction;
        direction = nextDirection;

        // Verify the true residual before acceptance and periodically to catch
        // loss of accuracy in the short recurrence. Restart only when needed.
        if(residualEstimate <= metricTolerance || nextBeta == 0 ||
           (iteration+1)%10 == 0 || iteration+1 == maxIterations) {
          State residual(b);
          if(recomputeResidual) residual = recomputeResidual(x);
          else residual -= Q(x);
          const double trueNorm = checkedNorm(residual);
          if(trueNorm <= tolerance) {
            x.toReal();
            logging::entry() << "MINRES ended after " << iteration+1
                             << " iterations; true relative residual=" << trueNorm/scale
                             << "; residual restarts=" << restarts << std::endl;
            return x;
          }
          State weighted = applyPreconditioner(residual);
          const double trueMetricNorm = metricNorm(residual, weighted);
          if(nextBeta == 0 || trueMetricNorm > 1.1*residualEstimate) {
            z = residual;
            z /= trueMetricNorm;
            v = weighted;
            v /= trueMetricNorm;
            previousZ *= 0;
            direction *= 0;
            previousDirection *= 0;
            beta = dbar = epsilon = sine = 0;
            cosine = -1;
            residualEstimate = trueMetricNorm;
            ++restarts;
            continue;
          }
        }
        previousZ = z;
        z = nextZ;
        z /= nextBeta;
        v = nextV;
        v /= nextBeta;
        beta = nextBeta;
        // if((iteration+1)%50 == 0)
          logging::entry() << "MINRES iteration " << iteration+1
                           << (precondition ? " preconditioned residual=" : " residual=")
                           << residualEstimate/metricScale << std::endl;
      }
      throw std::runtime_error("MINRES did not satisfy the requested tolerance");
    }
  }
}

#endif
