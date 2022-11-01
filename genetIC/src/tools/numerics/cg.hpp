#ifndef IC_CG_HPP
#define IC_CG_HPP

#include <functional>
#include <src/simulation/field/field.hpp>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/logging.hpp>

namespace tools {
  namespace numerics {

    //! Solve linear equation Qx = b, and return x, using conjugate gradient
    template<typename T>
    fields::Field<T> conjugateGradient(std::function<fields::Field<T>(const fields::Field<T> &)> Q,
                                       const fields::Field<T> &b,
                                       double rtol = 1e-6,
                                       double atol = 1e-12) {
      fields::Field<T> residual(b);
      fields::Field<T> direction = -residual;
      fields::Field<T> x = fields::Field<T>(b.getGrid(), false);

      double scale = b.norm();

      if(scale==0.0) {
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
        if (norm < rtol * scale || norm < atol)
          break;

        logging::entry() << "Conjugate gradient iteration " << i << " residual=" << norm << std::endl;

        // update direction for next cycle; must be Q-orthogonal to all previous updates
        double beta = residual.innerProduct(Q_direction) / direction.innerProduct(Q_direction);
        direction*=beta;
        direction-=residual;

      }
      logging::entry() << "Conjugate gradient ended after " << i << " iterations" << std::endl;

      return x;

    }

    template<typename T>
    fields::OutputField<T> conjugateGradient2(
        std::function<fields::OutputField<T>(fields::OutputField<T> &)> Q,
        fields::OutputField<T> &b,
        double rtol = 1e-6,
        double atol = 1e-12
    ) {
      fields::OutputField<T> residual(b);
      fields::OutputField<T> direction(b);
      fields::OutputField<T> x(b);
      
      double scale = 0;
      size_t dimension = 0;
      size_t Nlevel = b.getNumLevels();

      direction *= -1;
      x *= 0;
  
      for(auto ilevel = 0; ilevel<Nlevel; ++ilevel){
        auto field = b.getFieldForLevel(ilevel);
        assert(field.isFourier() == false);
        scale += field.norm();
        dimension += field.getGrid().size3;
      }

      if(scale==0.0) {
        logging::entry(logging::warning) << "Conjugate gradient: result is zero!" << std::endl;
        return x;
      }

      size_t iter = 0;
      for(; iter<dimension+1; ++iter) {
        auto Q_direction = Q(direction);
        double alpha = 0;
        for (auto ilevel = 0; ilevel<Nlevel; ++ilevel){
          // distance to travel in specified direction
          auto & resField = residual.getFieldForLevel(ilevel);
          auto & dirField = direction.getFieldForLevel(ilevel);
          auto & Q_dirField = Q_direction.getFieldForLevel(ilevel);
          assert(resField.isFourier() == false);
          assert(dirField.isFourier() == false);
          assert(Q_dirField.isFourier() == false);
          alpha += -resField.innerProduct(dirField) / dirField.innerProduct(Q_dirField);
        }

        x.addScaled(direction, alpha);

        residual = Q(x);

        residual.addScaled(b, -1);
        residual.toReal();
        double norm = 0;
        for (auto ilevel = 0; ilevel<Nlevel; ++ilevel){
          norm += residual.getFieldForLevel(ilevel).norm();
        }

        if (norm < rtol * scale || norm < atol)
          break;

        logging::entry() << "Conjugate gradient iteration " << iter << " residual=" << norm << std::endl;

        // update direction for next cycle; must be Q-orthogonal to all previous updates
        double beta = 0;
        for (auto ilevel = 0; ilevel<Nlevel; ++ilevel){
          auto & resField = residual.getFieldForLevel(ilevel);
          auto & dirField = direction.getFieldForLevel(ilevel);
          auto & Q_dirField = Q_direction.getFieldForLevel(ilevel);
          assert(resField.isFourier() == false);
          assert(dirField.isFourier() == false);
          assert(Q_dirField.isFourier() == false);
          beta += resField.innerProduct(Q_dirField) / dirField.innerProduct(Q_dirField);
        }

        direction *= beta;
        direction.addScaled(residual, -1);
      }
      logging::entry() << "Conjugate gradient ended after " << iter << " iterations" << std::endl;

      return x;
    }


    /* Adapted from https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf */
    template<typename T>
    fields::OutputField<T> minres(
        std::function<fields::OutputField<T>(fields::OutputField<T> &)> A,
        fields::OutputField<T> &b,
        double rtol = 1e-6,
        double atol = 1e-12
    ) {
      fields::OutputField<T> x(b);
      x *= 0;
      fields::OutputField<T> r(b);
      fields::OutputField<T> s(b);
      s = A(r);

      fields::OutputField<T> p(r);
      fields::OutputField<T> q(s);

      auto innerProduct = [](fields::OutputField<T> &a, fields::OutputField<T> &b) -> double {
        double result = 0;
        for (auto ilevel = 0; ilevel < a.getNumLevels(); ++ilevel) {
          result += a.getFieldForLevel(ilevel).innerProduct(b.getFieldForLevel(ilevel));
        }
        return result;
      };
      double scale = sqrt(innerProduct(b, b));
      double rho = innerProduct(r, s);

      size_t dimension = 0;
      for(auto ilevel = 0; ilevel<b.getNumLevels(); ++ilevel){
        dimension += b.getFieldForLevel(ilevel).getGrid().size;
      }

      // Start iteration
      size_t iter = 0;
      for(; iter<dimension; ++iter) {
        // We have q = A(p), but no need to compute it again
        double alpha = rho / innerProduct(q, q);
        x.addScaled(p, alpha);
        r.addScaled(q, -alpha);

        double norm = sqrt(innerProduct(r, r));
        if (norm < rtol * scale || norm < atol)
          break;

        logging::entry() << "MINRES iteration " << iter << " residual=" << norm << std::endl;
        s = A(r);
        double rhobar = rho;
        rho = innerProduct(r, s);
        double beta = rho / rhobar;
        p *= beta;
        p += r;

        q *= beta;
        q += s;
      }

      logging::entry() << "MINRES ended after " << iter << " iterations" << std::endl;

      return x;
    }
    template<typename T>
    fields::Field<T> minresField(
        std::function<fields::Field<T>(fields::Field<T> &)> A,
        fields::Field<T> &b,
        double rtol = 1e-6,
        double atol = 1e-12
    ) {
      fields::Field<T> x(b);
      x *= 0;
      fields::Field<T> r(b);
      fields::Field<T> s(b);
      s = A(r);

      fields::Field<T> p(r);
      fields::Field<T> q(s);

      double scale = b.norm();
      double rho = r.innerProduct(s);

      size_t dimension = b.getGrid().size3;

      // Start iteration
      size_t iter = 0;
      for(; iter<dimension; ++iter) {
        // We have q = A(p), but no need to compute it again
        double alpha = rho / q.innerProduct(q);
        x.addScaled(p, alpha);
        r.addScaled(q, -alpha);

        double norm = r.norm();
        if (norm < rtol * scale || norm < atol)
          break;

        logging::entry() << "MINRES iteration " << iter << " residual=" << norm << std::endl;
        s = A(r);
        double rhobar = rho;
        rho = r.innerProduct(s);
        double beta = rho / rhobar;
        p *= beta;
        p += r;

        q *= beta;
        q += s;
      }

      logging::entry() << "MINRES ended after " << iter << " iterations" << std::endl;

      return x;
    }


    template<typename T>
    fields::OutputField<T> bicgstab(
        std::function<fields::OutputField<T>(fields::OutputField<T> &)> A,
        fields::OutputField<T> &b,
        double rtol = 1e-6,
        double atol = 1e-12,
        const fields::OutputField<T> *x0 = nullptr
    ) {
      fields::OutputField<T> x(b);
      fields::OutputField<T> r(b);
      fields::OutputField<T> rOld(b);
      fields::OutputField<T> rhat(r);

      double rho;
      double omega;
      double rhoOld = 1;
      double alpha = 1;
      double omegaOld = 1;

      x *= 0;
      fields::OutputField<T> v(x);
      fields::OutputField<T> p(x);
      if (x0 != nullptr) {
        x = *x0;
        r.addScaled(A(x), -1);
      }

      auto innerProduct = [](fields::OutputField<T> &a, fields::OutputField<T> &b) -> double {
        double result = 0;
        for (auto ilevel = 0; ilevel < a.getNumLevels(); ++ilevel) {
          result += a.getFieldForLevel(ilevel).innerProduct(b.getFieldForLevel(ilevel));
        }
        return result;
      };

      assert (innerProduct(rhat, r) != 0);

      size_t dimension = 0;
      for(auto ilevel = 0; ilevel<b.getNumLevels(); ++ilevel){
        dimension += b.getFieldForLevel(ilevel).getGrid().size;
      }

      double scale = sqrt(innerProduct(b, b));

      // Start iteration
      size_t iter = 0;
      for(; iter<dimension; ++iter) {
        // rho_i = rhat_0 . r_{i-1}
        rho = innerProduct(rhat, rOld);
        // beta = (rho_i / rho_{i-1}) (alpha / omega_{i-1})
        double beta = (rho / rhoOld) * (alpha / omegaOld);
        // p_i = r_{i-1} + beta * (p_{i-1} - omega_{i-1} * v_{i-1})
        p *= beta;
        p.addScaled(rOld, 1);
        p.addScaled(v, -omegaOld * beta);
        // v_i = A.p_i
        v = A(p);
        // alpha = rho_i / (rhat_0 . v_i)
        alpha = rho / innerProduct(rhat, v);
        // x_{i} = x_{i-1} + alpha p_i
        x.addScaled(p, alpha);
        // s = r_{i-1} - alpha v_i
        auto s = rOld;
        s.addScaled(v, -alpha);
        double norm = std::sqrt(innerProduct(s, s));
        if (norm < rtol * scale || norm < atol) {
          break;
        }
        // t = A.s
        auto t = A(s);
        // omega_i = (t . s) / (t . t)
        omega = innerProduct(t, s) / innerProduct(t, t);
        // x_i = x_i + omega_i s
        x.addScaled(s, omega);
        // r_i = s - omega_i t
        r = s;
        r.addScaled(t, -omega);
        norm = std::sqrt(innerProduct(r, r));
        if (norm < rtol * scale || norm < atol)
          break;

        logging::entry() << "BiCGSTAB iteration " << iter << " residual=" << norm << std::endl;
        // Save old values
        rOld = r;
        rhoOld = rho;
        omegaOld = omega;
      }
      logging::entry() << "BiCGSTAB ended after " << iter << " iterations" << std::endl;
      return x;
    }
  }
}

#endif