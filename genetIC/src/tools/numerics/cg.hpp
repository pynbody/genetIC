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
     */
    template<typename T>
    fields::OutputField<T> conjugateGradient(
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> Q,
        const fields::OutputField<T> &b,
        double rtol = 1e-6, double atol = 1e-12,
        std::function<fields::OutputField<T>(const fields::OutputField<T> &)> recomputeResidual = {},
        size_t maxIterations = 0) {
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
        logging::entry() << "Conjugate gradient ended after 0 iterations" << std::endl;
        return x;
      }
      if(!maxIterations) {
        size_t dimension = 0;
        for(size_t level=0; level<b.getNumLevels(); ++level)
          dimension += b.getContext().getGridForLevel(level).size3;
        maxIterations = dimension+1;
      }
      for(size_t iteration=0; iteration<maxIterations; ++iteration) {
        const auto image = Q(direction);
        const double curvature = innerProduct(direction, image);
        if(!(curvature > 0) || !std::isfinite(curvature))
          throw std::runtime_error("Conjugate gradient breakdown before satisfying tolerance");
        const double step = rr/curvature;
        x.addScaled(direction, step);
        residual.addScaled(image, -step);
        double next = innerProduct(residual, residual);
        if(!std::isfinite(next) || next < 0)
          throw std::runtime_error("Conjugate gradient residual became nonfinite");
        // One operator application per ordinary iteration. Recompute the true
        // residual periodically and before acceptance; restart after replacement.
        if(std::sqrt(next) <= tolerance || (iteration+1)%50 == 0) {
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
                             << " iterations; true relative residual=" << std::sqrt(next)/scale << std::endl;
            return x;
          }
          direction = State(residual);
        } else {
          direction *= next/rr;
          direction += residual;
        }
        rr = next;
        if((iteration+1)%50 == 0)
          logging::entry() << "Conjugate gradient iteration " << iteration+1
                           << " residual=" << std::sqrt(rr)/scale << std::endl;
      }
      throw std::runtime_error("Conjugate gradient did not satisfy the requested tolerance");
    }
  }
}

#endif
