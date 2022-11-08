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
          const auto & left = a.getFieldForLevel(ilevel);
          const auto & right = b.getFieldForLevel(ilevel);
          result += left.innerProduct(right);
        }
        return result;
      };

      auto norm = [&innerProduct](fields::OutputField<T> &a) -> double {
        return std::sqrt(innerProduct(a, a));
      };

      double scale = sqrt(innerProduct(b, b));
      double rho = innerProduct(r, s);

      size_t dimension = 0;
      for(auto ilevel = 0; ilevel<b.getNumLevels(); ++ilevel){
        dimension += b.getFieldForLevel(ilevel).getGrid().size;
      }

      dimension *= 10;
      double old_norm = 0;
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

        // if (std::abs(norm - old_norm) / norm < rtol) {
        //   logging::entry(logging::warning) << "MINRES: stagnation detected" << std::endl;
        //   break;
        // }

        s = A(r);
        double rhobar = rho;
        rho = innerProduct(r, s);
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

        /* Adapted from https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf */
    template<typename T>
    fields::OutputField<T> minresYolo(
        std::function<fields::OutputField<T>(fields::OutputField<T> &)> A,
        fields::OutputField<T> &b,
        double rtol = 1e-6,
        double atol = 1e-12,
        fields::OutputField<T>* x0 = nullptr
    ) {
      assert (b.isRealOnAllLevels());

      int Nlevel = b.getNumLevels();

      auto innerProduct = [](fields::OutputField<T> &a, fields::OutputField<T> &b) -> double {
        double result = 0;
        for (auto ilevel = 0; ilevel < a.getNumLevels(); ++ilevel) {
          const auto & left = a.getFieldForLevel(ilevel);
          const auto & right = b.getFieldForLevel(ilevel);
          result += left.innerProduct(right);
        }
        return result;
      };

      auto norm = [&innerProduct](fields::OutputField<T> &a) -> double {
        return std::sqrt(innerProduct(a, a));
      };

      // my_minres(A, b, maxiter=None, x0=None, tol=1e-5):
      fields::OutputField<T> x(b.getContext(), b.getTransferType());
      x.getFieldForLevel(0);  // Trigger generation of grid
      x.toReal();
      if (x0 != nullptr) 
          x += *x0;

      int n = 0;

      for(auto ilevel = 0; ilevel<Nlevel; ++ilevel){
        const auto & field = b.getFieldForLevel(ilevel);
        assert(field.isFourier() == false);
        n += field.getGrid().size3;
      }

      auto maxiter = 5 * n;

      int istop = 0;
      int itn = 0;
      double Anorm = 0;
      double Acond = 0;
      double rnorm = 0;
      double ynorm = 0;

      double eps = std::numeric_limits<double>::epsilon();

      // Set up y and v for the first Lanczos vector v1.
      // y  =  beta1 P' v1,  where  P = C**(-1).
      // v is really P' v1.
      fields::OutputField<T> r1(b);
      if (x0 != nullptr)
          r1 -= A(*x0);
      fields::OutputField<T> y(r1);

      double beta1 = innerProduct(r1, y);

      if (beta1 < 0)
          throw std::runtime_error("minres: indefinite preconditioner");
      else if (beta1 == 0)
          return x;

      if (norm(b) == 0) {
          return b;
      }

      beta1 = std::sqrt(beta1);

      // Initialize other quantities
      double oldb = 0;
      double beta = beta1;
      double dbar = 0;
      double epsln = 0;
      double qrnorm = beta1;
      double phibar = beta1;
      double rhs1 = beta1;
      double rhs2 = 0;
      double tnorm2 = 0;
      double gmax = 0;
      double gmin = std::numeric_limits<double>::max();
      double cs = -1;
      double sn = 0;
      fields::OutputField<T> w(b.getContext(), b.getTransferType());
      fields::OutputField<T> w2(b.getContext(), b.getTransferType());
      w.getFieldForLevel(0);  // Trigger generation of grid
      w2.getFieldForLevel(0);  // Trigger generation of grid
      w.toReal();
      w2.toReal();
      fields::OutputField<T> r2(r1);
      for (auto itn = 1; itn <= maxiter; ++itn) {
          double s = 1.0/beta;
          // std::cout << "s = " << s << std::endl;
          auto v = y;
          v *= s;
          // std::cout << "v[0] = " << v.getFieldForLevel(0).getDataVector()[0] << std::endl;

          y = A(v); // TODO: spare some memory by storing v in y
          // std::cout << "y[0] = " << y.getFieldForLevel(0).getDataVector()[0] << std::endl;

          if (itn >= 2)
              y.addScaled(r1, -beta/oldb);

          double alfa = innerProduct(v, y);
          // std::cout << "alfa = " << alfa << std::endl;
          y.addScaled(r2, -alfa/beta);
          // std::cout << "y[0] = " << y.getFieldForLevel(0).getDataVector()[0] << std::endl;
          r1 = r2;
          // std::cout << "r1[0] = " << r1.getFieldForLevel(0).getDataVector()[0] << std::endl;
          r2 = y;
          // std::cout << "r2[0] = " << r2.getFieldForLevel(0).getDataVector()[0] << std::endl;
          y = r2;
          // std::cout << "y[0] = " << y.getFieldForLevel(0).getDataVector()[0] << std::endl;
          oldb = beta;
          // std::cout << "oldb = " << oldb << std::endl;
          beta = innerProduct(r2,y);
          // std::cout << "beta = " << beta << std::endl;
          if (beta < 0)
              throw std::runtime_error("non-symmetric matrix");
          beta = std::sqrt(beta);
          // std::cout << "beta = " << beta << std::endl;
          tnorm2 += alfa*alfa + oldb*oldb + beta*beta;
          // std::cout << "tnorm2 = " << tnorm2 << std::endl;

          if (itn == 1) {
              if (beta/beta1 <= 10*eps)
                  istop = -1;  // Terminate later
          }

          // Apply previous rotation Qk-1 to get
          //   [deltak epslnk+1] = [cs  sn][dbark    0   ]
          //   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

          double oldeps = epsln;
          double delta = cs * dbar + sn * alfa;   // delta1 = 0         deltak
          double gbar = sn * dbar - cs * alfa;   // gbar 1 = alfa1     gbar k
          epsln = sn * beta;     // epsln2 = 0         epslnk+1
          dbar = - cs * beta;   // dbar 2 = beta2     dbar k+1
          double root = std::sqrt(gbar*gbar + dbar*dbar);

          // Compute the next plane rotation Qk

          double gamma = std::sqrt(gbar*gbar + beta*beta);       // gammak
          gamma = std::max(gamma, eps);
          cs = gbar / gamma;             // ck
          sn = beta / gamma;             // sk
          double phi;
          phi = cs * phibar;             // phik
          phibar = sn * phibar;          // phibark+1

          // Update  x.
          double denom = 1.0/gamma;
          auto w1 = w2;
          w2 *= 0;
          w2 += w;
          w *= 0;
          w.addScaled(v, denom);
          w.addScaled(w1, -oldeps * denom);
          w.addScaled(w2, -delta * denom);

          // std::cout << "w[0] = " << w.getFieldForLevel(0).getDataVector()[0] << std::endl;

          x.addScaled(w, phi);
          // std::cout << "x[0] = " << x.getFieldForLevel(0).getDataVector()[0] << std::endl;

          // Go round again.

          gmax = std::max(gmax, gamma);
          gmin = std::min(gmin, gamma);
          double z = rhs1 / gamma;
          rhs1 = rhs2 - delta*z;
          rhs2 = - epsln*z;

          // Estimate various norms and test for convergence.

          Anorm = std::sqrt(tnorm2);
          ynorm = norm(x);
          double epsa = Anorm * eps;
          double epsx = Anorm * ynorm * eps;
          double epsr = Anorm * ynorm * rtol;
          double diag = gbar;

          if (diag == 0) diag = epsa;

          qrnorm = phibar;
          rnorm = qrnorm;
          double test1, test2;
          if (ynorm == 0 || Anorm == 0)
              test1 = 1e99;
          else
              test1 = rnorm / (Anorm*ynorm);    // ||r||  / (||A|| ||x||)
          if (Anorm == 0)
              test2 = 1e99;
          else
              test2 = root / Anorm;            // ||Ar|| / (||A|| ||r||)

          // Estimate  cond(A).
          // In this version we look at the diagonals of  R  in the
          // factorization of the lower Hessenberg matrix,  Q @ H = R,
          // where H is the tridiagonal matrix from Lanczos with one
          // extra row, beta(k+1) e_k^T.

          Acond = gmax/gmin;

          // See if any of the stopping criteria are satisfied.
          // In rare cases, istop is already -1 from above (Abar = const*I).

          if (istop == 0) {
              double t1 = 1 + test1;      // These tests work if tol < eps
              double t2 = 1 + test2;
              if (t2 <= 1)
                  istop = 2;
              if (t1 <= 1)
                  istop = 1;

              if (itn >= maxiter)
                  istop = 6;
              if (Acond >= 0.1/eps)
                  istop = 4;
              if (epsx >= beta1)
                  istop = 3;
              // if rnorm <= epsx   : istop = 2
              // if rnorm <= epsr   : istop = 1
              if (test2 <= rtol)
                  istop = 2;
              if (test1 <= rtol)
                  istop = 1;
          }

          // See if it is time to print something.
          bool prnt = (n <= 10) || (itn <= 10) || (itn >= maxiter-10) || ((itn % 10) == 0) || (qrnorm <= 10*epsx) || (qrnorm <= 10*epsr) || (Acond <= 1e-2/eps) || (istop != 0);

          if (prnt) {
              // str1 = '%6g %12.5e %10.3e' % (itn, x[0], test1)
              // str2 = ' %10.3e' % (test2,)
              // str3 = ' %8.1e %8.1e %8.1e' % (Anorm, Acond, gbar/Anorm)
            char buffer[100];
            sprintf(buffer, "%6d %12.5e %10.3e %10.3e %8.1e %8.1e %8.1e", itn, x.getFieldForLevel(0).getDataVector() [0], test1, test2, Anorm, Acond, gbar/Anorm);
            logging::entry() << "MINRES " << buffer << std::endl;

          }

          if (istop != 0)
              break;  // TODO check this
      }

      int info;
      if (istop == 6)
          info = maxiter;
      else
          info = 0;

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

      dimension *= 10;

      double scale = sqrt(innerProduct(b, b));
      double old_norm = 0;

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

        if (std::abs(norm - old_norm) / norm < rtol) {
          logging::entry(logging::warning) << "BiCGSTAB: stagnation detected" << std::endl;
          break;
        }
        old_norm = norm;

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