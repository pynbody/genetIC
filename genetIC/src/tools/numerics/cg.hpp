#ifndef IC_CG_HPP
#define IC_CG_HPP

#include <functional>
#include <src/simulation/field/field.hpp>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/logging.hpp>

namespace tools {
  namespace numerics {

    //! Solve linear equation Qx = b, and return x, using conjugate gradient
    template <typename T>
    fields::Field<T> conjugateGradient(std::function<fields::Field<T>(const fields::Field<T> &)> Q,
                                       const fields::Field<T> &b,
                                       double rtol = 1e-4, // changed from 1e-6
                                       double atol = 1e-12)
    {
      fields::Field<T> residual(b);
      fields::Field<T> direction = -residual;
      fields::Field<T> x = fields::Field<T>(b.getGrid(), false);

      double scaleNorm = b.norm();
      double scaleMax = b.Maximum(direction);

      if(scaleNorm==0.0) {
        logging::entry(logging::warning) << "Conjugate gradient: result is zero!" << std::endl;
        return x;
      }

      size_t dimension = b.getGrid().size3;

      size_t i;

      for(i=0; i<dimension+1; ++i) {

        fields::Field<T> Q_direction = Q(direction);
        // distance to travel in specified direction
        double alpha = -residual.innerProduct(direction) / direction.innerProduct(Q_direction);
        x.addScaled(direction, alpha);

        residual = Q(x);
        residual -= b;

        auto norm = residual.norm();
        auto max = residual.Maximum(direction);
        // if (norm < rtol * scaleNorm || norm < atol)
        //   break;
        if (max < 1e-6 * scaleMax)
          break;

        logging::entry() << "Conjugate gradient iteration " << i << std::endl;
        logging::entry() << "conditionNorm=" << rtol * scaleNorm << " residual=" << norm << std::endl;
        logging::entry() << "conditionMax=" << rtol * scaleMax << " maximum=" << max << std::endl;

        // update direction for next cycle; must be Q-orthogonal to all previous updates
        double beta = residual.innerProduct(Q_direction) / direction.innerProduct(Q_direction);
        direction*=beta;
        direction-=residual;

      }
      logging::entry() << "Conjugate gradient ended after " << i << " iterations" << std::endl;

      return x;
    }

        /* Adapted from https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf */
    template<typename T>
    fields::Field<T> minres(
        std::function<fields::Field<T>(fields::Field<T> &)> A,
        fields::Field<T> &b,
        double rtol = 1e-4, // changed from 1e-6
        double atol = 1e-12)
    {
      fields::Field<T> x(b);
      x *= 0;
      fields::Field<T> r(b);
      fields::Field<T> s(b);
      s = A(r);

      fields::Field<T> p(r);
      fields::Field<T> q(s);

      auto norm = b.norm();

      double scale = sqrt(b.innerProduct(b));
      double rho = r.innerProduct(s);

      size_t dimension = b.getGrid().size3;

      dimension *= 10;
      double old_norm = 0;
      // Start iteration
      size_t iter = 0;
      for(; iter<dimension; ++iter) {
        // We have q = A(p), but no need to compute it again
        double alpha = rho / q.innerProduct(q);
        x.addScaled(p, alpha);
        r.addScaled(q, -alpha);

        double norm = sqrt(r.innerProduct(r));
        double max = r.Maximum(r);
        if (max < 1e-5 * scale || norm < atol)
          break;
        logging::entry() << "MINRES iteration " << iter << " residual=" << norm << std::endl;

        // if (std::abs(norm - old_norm) / norm < rtol) {
        //   logging::entry(logging::warning) << "MINRES: stagnation detected" << std::endl;
        //   break;
        // }

        s = A(r);
        double rhobar = rho;
        rho = r.innerProduct(s);
        double beta = rho / rhobar;
        p *= beta;
        p += r;

        q *= beta;
        q += s;

        // Save old norm
        old_norm = norm;
      }

      logging::entry() << "MINRES ended after " << iter << " iterations" << std::endl;

      return x;

    }
  }
}

#endif