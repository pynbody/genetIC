#ifndef _INTERPOLATION_HPP_INCLUDED
#define _INTERPOLATION_HPP_INCLUDED

#include <type_traits>

/*!
    \namespace numerics
    \brief Provides numerical methods such as interpolation or fourier transforms.

 */
namespace numerics {

  /** Spline interpolation object
     *
     * Wraps the GSL spline library
     */
  template<typename T>
  class Interpolator : std::enable_shared_from_this<Interpolator<T>> {

    // At present, we only have GSL interpolation which _requires_ double arguments
    static_assert(std::is_same<T, double>::value, "Only support interpolation over doubles");

  private:
    gsl_interp_accel *acc;
    gsl_spline *spline;

  public:
    Interpolator() {
      acc = nullptr;
      spline = nullptr;
    }

    Interpolator(std::vector<T> &x, std::vector<T> &y) : Interpolator() {
      initialise(x, y);
    }

    void initialise(std::vector<T> &x, std::vector<T> &y) {
      assert(x.size() == y.size());
      deinitialise();
      acc = gsl_interp_accel_alloc();
      spline = gsl_spline_alloc(gsl_interp_cspline, y.size());
      gsl_spline_init(spline, x.data(), y.data(), x.size());
    }

    void deinitialise() {
      if (spline != nullptr)
        gsl_spline_free(spline);
      if (acc != nullptr)
        gsl_interp_accel_free(acc);
      spline = nullptr;
      acc = nullptr;
    }

    ~Interpolator() {
      deinitialise();
    }

    T operator()(T x) const {
      return gsl_spline_eval(spline, x,
                             nullptr); // not able to pass accelerator as final argument because it is not thread-safe
    }


  };
}

#endif

