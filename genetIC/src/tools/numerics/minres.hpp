#ifndef IC_MINRES_HPP
#define IC_MINRES_HPP

#include <functional>
#include <src/simulation/field/field.hpp>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/logging.hpp>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "cg.hpp"

namespace tools {
  namespace numerics {

    /* Adapted from https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf */
    template<typename T>
    fields::OutputField<T> minres(std::function<fields::OutputField<T>(const fields::OutputField<T> &)> A,
                            const fields::OutputField<T> &b,
                            double rtol = 1e-6,
                            double atol = 1e-12 )
    {
      fields::OutputField<T> x(b.getContext(), b.getTransferType());
      x.getFieldForLevel(0);  // Trigger allocation
      fields::OutputField<T> r(b);
      fields::OutputField<T> s(b);
      s = A(r);

      fields::OutputField<T> p(r);
      fields::OutputField<T> q(s);

      double scale = tools::numerics::norm(b);
      double rho = tools::numerics::innerProduct(r, s);

      size_t dimension = 0;
      for (auto ilevel = 0; ilevel < b.getNumLevels(); ++ilevel) {
        const auto ctxt = r.getContext();
        dimension += ctxt.getGridForLevel(ilevel).size3;
      }

      size_t max_iterations = dimension * 10;
      double old_norm = 0;

      size_t iter = 0;
      
      for(; iter<max_iterations; ++iter) {
        // We have q = A(p), but no need to compute it again
        double alpha = rho / tools::numerics::innerProduct(q, q);
        x.addScaled(p, alpha);
        r.addScaled(q, -alpha);

        double norm = tools::numerics::norm(r);

        if (norm < rtol * scale || norm < atol)
          break;

        logging::entry() << "MINRES iteration " << iter << " residual=" << norm/scale << std::endl;

        s = A(r);
        double rhobar = rho;
        rho = tools::numerics::innerProduct(r, s);
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