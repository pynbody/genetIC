#ifndef IC_CG_HPP
#define IC_CG_HPP

#include <functional>
#include <src/simulation/field/field.hpp>
#include <src/tools/data_types/complex.hpp>
#include <src/tools/logging.hpp>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

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
      double scaleMax = b.Maximum();

      if(scaleNorm==0.0) {
        logging::entry(logging::warning) << "Conjugate gradient: result is zero!" << std::endl;
        return x;
      }

      auto toltest = 1e-6;

      logging::entry() << "conditionNorm=" << toltest * scaleNorm << " conditionMax=" << toltest * scaleMax << std::endl;

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
        auto max = residual.Maximum();
        // if (norm < rtol * scaleNorm || norm < atol)
        //   break;

        if (max < toltest * scaleMax)
          break;

        logging::entry() << "Conjugate gradient iteration " << i << "residual=" << norm << " maximum=" << max << std::endl;

        // update direction for next cycle; must be Q-orthogonal to all previous updates
        double beta = residual.innerProduct(Q_direction) / direction.innerProduct(Q_direction);
        direction*=beta;
        direction-=residual;

      }
      logging::entry() << "Conjugate gradient ended after " << i << " iterations" << std::endl;

      return x;
    }


    /* Adapted from https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf */
    /*
    template<typename T>
    fields::Field<T> minres(std::function<fields::Field<T>(fields::Field<T> &)> A,
                            fields::Field<T> &b,
                            double rtol = 1e-4, // changed from 1e-6
                            double atol = 1e-12,
                            bool restart = false,
                            std::string output_path = "")
    {
      fields::Field<T> x(b);
      x *= 0;
      fields::Field<T> r(b);
      fields::Field<T> s(b);
      s = A(r);

      fields::Field<T> p(r);
      fields::Field<T> q(s);

      auto toltest = 1e-7;

      double scaleNorm = sqrt(b.innerProduct(b));
      double scaleMax = b.Maximum();
      double rho = r.innerProduct(s);

      size_t dimension = b.getGrid().size3;

      dimension *= 10;
      double old_norm = 0;

      logging::entry() << "conditionNorm=" << toltest * scaleNorm << " conditionMax=" << toltest * scaleMax << std::endl;


      const double brakeTime = 0.1;  // Desired brake time in hours
      const std::time_t brakeDuration = brakeTime * 3600; // Calculate the duration in seconds for the brake time
      const std::time_t startTime = std::time(nullptr); // Get the current time
      logging::entry() << "conditionNorm=" << toltest * scaleNorm << " conditionMax=" << toltest * scaleMax << std::endl;

      // Start iteration
      size_t iter = 0;

      //setting variables from last restart
      if (restart == true) {
          logging::entry() << "Loading restart variables" << std::endl;

          std::string variables = ".variables";
          std::string directory = output_path + variables;

          x.loadGridData(directory + "/x");
          r.loadGridData(directory + "/r");
          q.loadGridData(directory + "/q");
          p.loadGridData(directory + "/p");
        
          std::ifstream file;

          file.open (directory + "/rho.txt");
          file >> rho;
          file.close();

          file.open (directory + "/iter.txt");
          file >> iter;
          file.close();

      }

      for(; iter<dimension; ++iter) {

        // We have q = A(p), but no need to compute it again
        double alpha = rho / q.innerProduct(q);
        x.addScaled(p, alpha);
        r.addScaled(q, -alpha);

        double norm = sqrt(r.innerProduct(r));
        double max = r.Maximum();

        if (max < toltest * scaleMax || norm < atol)
          break;

        const std::time_t currentTime = std::time(nullptr); // Get the current time
        const std::time_t elapsedTime = currentTime - startTime; // Calculate the elapsed time
        if (elapsedTime >= brakeDuration) {
          logging::entry() << "Maximum time reached" << std::endl;
          break;
        }
        
        logging::entry() << "MINRES iteration " << iter << std::endl;
        logging::entry() << "conditionNorm=" << toltest * scaleNorm << " residual=" << norm << std::endl;
        logging::entry() << "conditionMax=" << toltest * scaleMax << " maximum=" << max << std::endl;

        // if (std::abs(norm - old_norm) / norm < toltest) {
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
        
        //save current minimization variables if time limit reached
        if (elapsedTime >= brakeDuration) {
          logging::entry() << "Maximum time reached" << std::endl;
          
          std::string variables = ".variables";
          std::string directory = output_path + variables;

          mkdir(directory.c_str(), 0777);

          //save fields
          x.dumpGridData(directory + "/x");
          r.dumpGridData(directory + "/r");
          q.dumpGridData(directory + "/q");
          p.dumpGridData(directory + "/p");

          //save variables
          std::ofstream file;

          file.open (directory + "/rho.txt");
          file << rho;
          file.close();

          file.open (directory + "/iter.txt");
          file << iter + 1;
          file.close();

          break;
        }

        // Save old norm
        old_norm = norm;
      }

      logging::entry() << "MINRES ended after " << iter << " iterations" << std::endl;

      return x;

    }
    */
  }
}

#endif