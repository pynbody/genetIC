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

namespace tools {
  namespace numerics {

    /* Adapted from https://stanford.edu/group/SOL/reports/SOL-2011-2R.pdf */
    template<typename T>
    fields::Field<T> minres(std::function<fields::Field<T>(fields::Field<T> &)> A,
                            fields::Field<T> &b,
                            double rtol = 2e-4, // changed from 1e-6
                            double atol = 1e-12,
                            bool restart = false,
                            bool stop = false,
                            double brakeTime = 0,
                            std::string output_path="")
    {
      fields::Field<T> x(b);
      x *= 0;
      fields::Field<T> r(b);
      fields::Field<T> s(b);
      s = A(r);

      fields::Field<T> p(r);
      fields::Field<T> q(s);

      // auto toltest = 1e-8;     // Default tolerance (meant to even include 1024^3 potential splicing)

      double scaleNorm = sqrt(b.innerProduct(b));
      double scaleMax = b.Maximum();
      double rho = r.innerProduct(s);

      size_t dimension = b.getGrid().size3;

      /*
      if (rtol != 2e-4) {

        // Lowering tolerance for lower resolution fields
        if (dimension == 512*512*512)
          toltest = 1e-7;
        if (dimension == 256*256*256) {
          toltest = 1e-6;
          brakeTime = 24; // Perhaps a time dependence instead of trying to achieve the
                          // exact tolerance is better
        }
      }
      */

      size_t max_iterations = dimension * 10;
      double old_norm = 0;

      size_t iter = 0;
      //setting variables from last restart
      /*
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
      */

      if (stop == true){

        const std::time_t brakeDuration = brakeTime * 3600; // Calculate the duration in seconds for the brake time
        const std::time_t startTime = std::time(nullptr);   // Get the current time at splicing start

        logging::entry() << "Splicing will stop after " << brakeTime << " hours, or when the maximum drops below " << rtol * scaleMax << std::endl;

        // Start iteration
        for(; iter<max_iterations; ++iter) {

          // We have q = A(p), but no need to compute it again
          double alpha = rho / q.innerProduct(q);
          x.addScaled(p, alpha);
          r.addScaled(q, -alpha);

          double norm = sqrt(r.innerProduct(r));
          double max = r.Maximum();

          if (max < rtol * scaleMax || norm < atol)
            break;

          const std::time_t currentTime = std::time(nullptr); // Get the current time
          const std::time_t elapsedTime = currentTime - startTime; // Calculate the elapsed time
          
          logging::entry() << "MINRES iteration " << iter << " maximum=" << max << std::endl;

          s = A(r);
          double rhobar = rho;
          rho = r.innerProduct(s);
          double beta = rho / rhobar;
          p *= beta;
          p += r;

          q *= beta;
          q += s;

          if (elapsedTime >= brakeDuration) {

            logging::entry() << "Maximum time reached" << std::endl;
            break;
          }
            //if (max < rtol * scaleMax || norm < atol) {

            //  logging::entry() << "Maximum time reached" << std::endl;
            
              // Save current minimization variables if time limit reached
              /*            
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
              */

            //  break;

            //}
        }

          // Save old norm
          // old_norm = norm;
      }
      else {

        logging::entry() << "Splicing will stop when the maximum drops below " << rtol * scaleMax << std::endl;
        for(; iter<max_iterations; ++iter) {

          // We have q = A(p), but no need to compute it again
          double alpha = rho / q.innerProduct(q);
          x.addScaled(p, alpha);
          r.addScaled(q, -alpha);

          double norm = sqrt(r.innerProduct(r));
          double max = r.Maximum();

          if (max < rtol * scaleMax || norm < atol)
            break;
          
          logging::entry() << "MINRES iteration " << iter << " maximum=" << max << std::endl;

          s = A(r);
          double rhobar = rho;
          rho = r.innerProduct(s);
          double beta = rho / rhobar;
          p *= beta;
          p += r;

          q *= beta;
          q += s;
        }
      }

      logging::entry() << "MINRES ended after " << iter << " iterations" << std::endl;

      return x;
    }
  }
}

#endif