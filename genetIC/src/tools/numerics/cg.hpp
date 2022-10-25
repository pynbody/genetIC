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
      double r2 = innerProduct(r, r);
      double rho = innerProduct(r, s);

      size_t dimension = 0;
      for(auto ilevel = 0; ilevel<b.getNumLevels(); ++ilevel){
        dimension += b.getFieldForLevel(ilevel).getGrid().size3;
      }

      // Start iteration
      size_t iter = 0;
      for(; iter<dimension; ++iter) {
        // We have q = A(p), but no need to compute it again
        double alpha = rho / innerProduct(q, q);
        x.addScaled(p, alpha);
        r.addScaled(q, -alpha);

        double norm = innerProduct(r, r);
        if (norm < rtol * r2 || norm < atol)
          break;

        logging::entry() << "Conjugate gradient iteration " << iter << " residual=" << norm << std::endl;
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
  }
}

#endif