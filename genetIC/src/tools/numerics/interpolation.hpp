#ifndef _INTERPOLATION_HPP_INCLUDED
#define _INTERPOLATION_HPP_INCLUDED

#include <type_traits>
#include <gsl/gsl_spline.h>

namespace tools {
/*!
    \namespace tools::numerics
    \brief Provides numerical methods such as interpolation or fourier transforms.

 */
  namespace numerics {

    /*! \class Interpolator
        \brief Spline interpolation object
       *
       * Wraps the GSL spline library
       */
    template<typename T>
    class Interpolator : std::enable_shared_from_this<Interpolator<T>> {

      // At present, we only have GSL interpolation which _requires_ double arguments
      static_assert(std::is_same<T, double>::value, "Only support interpolation over doubles");

    private:
      gsl_interp_accel *acc; //!< Accelerator which allows rapid searching to find the right polynomial to use
      gsl_spline *spline; //!< Spline object that performs the evaluation

    protected:
      T minx, maxx; //!< Minimum and maximum permitted values of function argument

    public:
      //! Default constructor - no spline or acc specified
      Interpolator() {
        acc = nullptr;
        spline = nullptr;
      }

      //! Constructor defined from a source and value vector
      Interpolator(std::vector<T> &x, std::vector<T> &y) : Interpolator() {
        initialise(x, y);
      }

      //! Initialise the gsl spline required for these data
      virtual void initialise(std::vector<T> &x, std::vector<T> &y) {
        assert(x.size() == y.size());
        deinitialise();
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, y.size());
        gsl_spline_init(spline, x.data(), y.data(), x.size());
        minx = x[0];
        maxx = x.back();
      }

      //! De-initialise the gsl spline:
      virtual void deinitialise() {
        if (spline != nullptr)
          gsl_spline_free(spline);
        if (acc != nullptr)
          gsl_interp_accel_free(acc);
        spline = nullptr;
        acc = nullptr;
      }

      //! Destructor (de-initialises the spline)
      ~Interpolator() {
        deinitialise();
      }

      //! Evaluates the spline-function at the specified value of x
      virtual T operator()(T x) const {
        if(x<minx || x>maxx) return T(0);
        return gsl_spline_eval(spline, x,
                               nullptr); // not able to pass accelerator as final argument because it is not thread-safe
      }


    };

    template<typename T>
    class LogInterpolator : Interpolator<T> {
    protected:
      std::vector<T> logx;
      std::vector<T> logy;

    public:
      void initialise(std::vector<T> &x, std::vector<T> &y) override {
        const T inf = std::numeric_limits<T>::infinity();
        assert(x.size()==y.size());
        for(size_t i = 0; i<x.size(); ++i) {
          T logx_i = std::log(x[i]);
          T logy_i = std::log(y[i]);
          if(logx_i!=inf && logx_i!=-inf && !std::isnan(logx_i)
            && logy_i!=inf && logy_i!=-inf && !std::isnan(logy_i)) {
            logx.push_back(logx_i);
            logy.push_back(logy_i);
          }
        }
        Interpolator<T>::initialise(logx, logy);
      }

      T operator()(T x) const override {
        T logx = std::log(x);
        if(logx<this->minx || logx>this->maxx) return T(0);
        return std::exp(Interpolator<T>::operator()(logx));
      }
    };
  }
}
#endif

