#ifndef IC_CG_HPP
#define IC_CG_HPP

#include <functional>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <src/simulation/field/multilevelfield.hpp>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/logging.hpp>

namespace tools {
  namespace numerics {

    template<typename T>
    T innerProduct(const fields::OutputField<T> &a, const fields::OutputField<T> &b) {
      T result = 0;
      for (size_t ilevel = 0; ilevel < a.getNumLevels(); ++ilevel) {
        const auto& left = a.getFieldForLevel(ilevel);
        const auto& right = b.getFieldForLevel(ilevel);
        result += left.innerProduct(right);
      }
      return result;
    }

    template<typename T>
    double norm(const fields::OutputField<T> &a) {
      return std::sqrt(innerProduct(a, a));
    }

    /*! Solve Qx=b by conjugate gradient. Q must be positive definite on the
     * reachable subspace (consistent positive-semidefinite systems also work).
     * recomputeResidual optionally evaluates b-Qx through the actual returned
     * observable, avoiding cancellation between separately evaluated operators.
     * It must use the same sign as b-Qx. Zero maxIterations uses dimension+1.
     * precondition applies an SPD approximate inverse of Q. Acceptance always
     * uses the unpreconditioned residual, in the original units.
     */
    template<typename T>
    fields::OutputField<T> conjugateGradient(
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> Q,
        const fields::OutputField<T> &b,
        double rtol = 1e-6, double atol = 1e-12,
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> recomputeResidual = {},
        size_t maxIterations = 0,
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> precondition = {}) {
      using State = fields::OutputField<T>;
      if(!std::isfinite(rtol) || !std::isfinite(atol) || rtol < 0 || atol < 0)
        throw std::runtime_error("Conjugate gradient tolerances must be finite and nonnegative");
      State residual(b), direction(b);
      State x(b.getContext(), b.getTransferType());
      x.getFieldForLevel(0);
      double rr = innerProduct(residual, residual);
      if(!std::isfinite(rr) || rr < 0)
        throw std::runtime_error("Conjugate gradient received a nonfinite residual");
      const double scale = std::sqrt(rr);
      const double tolerance = std::max(rtol*scale, atol);
      if(scale <= tolerance) {
        // Match the real-space result returned after an ordinary CG update.
        x.toReal();
        logging::entry() << "Conjugate gradient ended after 0 iterations" << std::endl;
        return x;
      }
      if(!maxIterations) {
        size_t dimension = 0;
        for(size_t level=0; level<b.getNumLevels(); ++level)
          dimension += b.getContext().getGridForLevel(level).size3;
        maxIterations = dimension+1;
      }
      State weighted = precondition ? precondition(residual) : State(residual);
      double rho = innerProduct(residual, weighted);
      direction = weighted;
      size_t restarts = 0;
      for(size_t iteration=0; iteration<maxIterations; ++iteration) {
        if(!(rho > 0) || !std::isfinite(rho))
          throw std::runtime_error("Conjugate gradient preconditioner must be positive definite");
        const auto image = Q(direction);
        const double curvature = innerProduct(direction, image);
        if(!(curvature > 0) || !std::isfinite(curvature))
          throw std::runtime_error("Conjugate gradient breakdown before satisfying tolerance");
        const double step = rho/curvature;
        x.addScaled(direction, step);
        residual.addScaled(image, -step);
        double next = innerProduct(residual, residual);
        if(!std::isfinite(next) || next < 0)
          throw std::runtime_error("Conjugate gradient residual became nonfinite");
        bool restart = false;
        // Check actual residuals periodically and before acceptance. Discarding
        // conjugacy every 50 iterations needlessly slows well-resolved solves;
        // restart only when the recurrence has lost significant accuracy.
        if(std::sqrt(next) <= tolerance || (iteration+1)%10 == 0 || iteration+1 == maxIterations) {
          State recursiveResidual(residual);
          if(recomputeResidual) residual = recomputeResidual(x);
          else {
            residual = State(b);
            residual -= Q(x);
          }
          next = innerProduct(residual, residual);
          if(!std::isfinite(next) || next < 0)
            throw std::runtime_error("Conjugate gradient true residual became nonfinite");
          if(std::sqrt(next) <= tolerance) {
            logging::entry() << "Conjugate gradient ended after " << iteration+1
                             << " iterations; true relative residual=" << std::sqrt(next)/scale
                             << "; residual restarts=" << restarts << std::endl;
            return x;
          }
          recursiveResidual -= residual;
          restart = norm(recursiveResidual) > 0.1*std::sqrt(next);
        }
        weighted = precondition ? precondition(residual) : State(residual);
        const double nextRho = innerProduct(residual, weighted);
        if(restart) {
          direction = weighted;
          ++restarts;
        } else {
          direction *= nextRho/rho;
          direction += weighted;
        }
        rr = next;
        rho = nextRho;
        // if((iteration+1)50 == 0)
          logging::entry() << "Conjugate gradient iteration " << iteration+1
                           << " residual=" << std::sqrt(rr)/scale << std::endl;
      }
      throw std::runtime_error("Conjugate gradient did not satisfy the requested tolerance");
    }
  }
}

#endif
