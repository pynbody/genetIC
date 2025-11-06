#ifndef IC_CG_HPP
#define IC_CG_HPP

#include <functional>
#include <src/simulation/field/field.hpp>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/logging.hpp>

namespace tools {
  namespace numerics {

    template<typename T>
    T innerProduct(const fields::OutputField<T> &a, const fields::OutputField<T> &b) {
      T result = 0;
      for (auto ilevel = 0; ilevel < a.getNumLevels(); ++ilevel) {
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

    //! Solve linear equation Qx = b, and return x, using conjugate gradient
    template<typename T>
    fields::OutputField<T> conjugateGradient(std::function<fields::OutputField<T>(const fields::OutputField<T> &)> Q,
                                       const fields::OutputField<T> &b,
                                       double rtol = 1e-6,
                                       double atol = 1e-12) {
      fields::OutputField<T> residual(b);
      fields::OutputField<T> direction(residual);
      direction *= -1;
      fields::OutputField<T> x = fields::OutputField<T>(b.getContext(), b.getTransferType());
      x.getFieldForLevel(0); // trigger allocation

      double scale = norm(residual);

      if(scale==0.0) {
        logging::entry(logging::warning) << "Conjugate gradient: result is zero!" << std::endl;
        return x;
      }

      size_t dimension = 0;
      for (auto ilevel = 0; ilevel < b.getNumLevels(); ++ilevel) {
        const auto ctxt = residual.getContext();
        dimension += ctxt.getGridForLevel(ilevel).size3;
      }

      size_t i;

      for(i=0; i<dimension+1; ++i) {

        auto Q_direction = Q(direction);
  
        // distance to travel in specified direction
        double alpha = -innerProduct(residual, direction) / innerProduct(direction, Q_direction);

        x.addScaled(direction, alpha);

        residual = Q(x);
        residual -= b;

        auto res_norm = norm(residual);
        if (res_norm < rtol * scale || res_norm < atol)
          break;

        logging::entry() << "Conjugate gradient iteration " << i << " residual=" << res_norm << "/" << scale << std::endl;

        // update direction for next cycle; must be Q-orthogonal to all previous updates
        double beta = innerProduct(residual, Q_direction) / innerProduct(direction, Q_direction);
        direction*=beta;
        direction-=residual;

      }
      logging::entry() << "Conjugate gradient ended after " << i << " iterations" << std::endl;

      return x;

    }
  }
}

#endif