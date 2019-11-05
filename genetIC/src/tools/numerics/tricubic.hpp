//
// Created by Andrew Pontzen on 2019-09-15.
//

#ifndef IC_TRICUBIC_HPP
#define IC_TRICUBIC_HPP

namespace numerics {
  /* \brief Fast calculation of a^n for positive integer n
   * 
   * NB this routine is *much* faster than calling std::pow, but ONLY if the compiler optimizes using the
   * constexpr nature of the routine. This seems to happen on recent versions of Clang and GCC at -O2 or above.
   */
  template<typename T>
  constexpr inline T fastpow(T a, unsigned int n) {
    T result = 1.0;
    for(unsigned int i=0; i<n; i++) {
      result*=a;
    }
    return result;
  }

  // If complex interpolation is attempted, we need to be able to multiply complex<T> by an integer
  template<typename T>
  std::complex<T> operator*(int a, std::complex<T> b) {
    return std::complex<T>(b.real()*a, b.imag()*a);
  }

  /* \brief A class to perform local function estimation by tricubic interpolation on the domain [0,1]^3
   *
   * That is, the function is approximated by \sum_ijk=0^3 a_{ijk} x^i y^j z^k
   *
   * The 64 coefficients a_{ijk} are conceptually found by specifying the values, all first derivatives, plus
   * d2/dxdy, d2/dxdz, d2/dydz, d3/dxdydz at each corner of the unit cube. However in practice these derivatives
   * are estimated numerically from the values of the function at the corners of the cube [-1,2]^3. This generates
   * 4^3=64 constraints as expected.
   *
   */
  template<typename T>
  class LocalUnitTricubicApproximation {
  protected:
    T a000, a001, a002, a003, a010, a011, a012, a013, a020, a021, a022, \
      a023, a030, a031, a032, a033, a100, a101, a102, a103, a110, a111, \
      a112, a113, a120, a121, a122, a123, a130, a131, a132, a133, a200, \
      a201, a202, a203, a210, a211, a212, a213, a220, a221, a222, a223, \
      a230, a231, a232, a233, a300, a301, a302, a303, a310, a311, a312, \
      a313, a320, a321, a322, a323, a330, a331, a332, a333;
  public:
    /* Initialise using samples from the function.
     *
     * The indices are such that cellValues[1][1][1] gives the lower-left corner of the unit cube (0,1)^3,
     * cellValues[2][2][2] gives its upper right corner. Data outside the unit cube (i.e. with indices 0 or 3)
     * is used to numerically estimated derivatives. */
    explicit LocalUnitTricubicApproximation(T cellValues[4][4][4]) {
      initCoeffsFromCellValues(cellValues);
    }

    LocalUnitTricubicApproximation() {}

    // Evaluate the interpolated function within the unit cube [0,1]^3
    T operator()(T x, T y, T z) const {
      assert(x>=0 && x<1 && y>=0 && y<1 && z>=0 && z<1);
      // If speed-ups are necessary here, the values of x,y,z for which this is called are likely to be highly
      // predictable and results globally cachable. Caching the a... coefficients could also
      // be undertaken (but not here, obviously).
      return a000 + a100*x + a200*fastpow(x,2) + a300*fastpow(x,3) + a010*y + a110*x*y + a210*fastpow(x,2)*y + a310*fastpow(x,3)*y + a020*fastpow(y,2) + a120*x*fastpow(y,2) +
             a220*fastpow(x,2)*fastpow(y,2) + a320*fastpow(x,3)*fastpow(y,2) + a030*fastpow(y,3) + a130*x*fastpow(y,3) + a230*fastpow(x,2)*fastpow(y,3) + a330*fastpow(x,3)*fastpow(y,3) + a001*z +
             a101*x*z + a201*fastpow(x,2)*z + a301*fastpow(x,3)*z + a011*y*z + a111*x*y*z + a211*fastpow(x,2)*y*z + a311*fastpow(x,3)*y*z + a021*fastpow(y,2)*z + a121*x*fastpow(y,2)*z +
             a221*fastpow(x,2)*fastpow(y,2)*z + a321*fastpow(x,3)*fastpow(y,2)*z + a031*fastpow(y,3)*z + a131*x*fastpow(y,3)*z + a231*fastpow(x,2)*fastpow(y,3)*z + a331*fastpow(x,3)*fastpow(y,3)*z +
             a002*fastpow(z,2) + a102*x*fastpow(z,2) + a202*fastpow(x,2)*fastpow(z,2) + a302*fastpow(x,3)*fastpow(z,2) + a012*y*fastpow(z,2) + a112*x*y*fastpow(z,2) +
             a212*fastpow(x,2)*y*fastpow(z,2) + a312*fastpow(x,3)*y*fastpow(z,2) + a022*fastpow(y,2)*fastpow(z,2) + a122*x*fastpow(y,2)*fastpow(z,2) + a222*fastpow(x,2)*fastpow(y,2)*fastpow(z,2) +
             a322*fastpow(x,3)*fastpow(y,2)*fastpow(z,2) + a032*fastpow(y,3)*fastpow(z,2) + a132*x*fastpow(y,3)*fastpow(z,2) + a232*fastpow(x,2)*fastpow(y,3)*fastpow(z,2) +
             a332*fastpow(x,3)*fastpow(y,3)*fastpow(z,2) + a003*fastpow(z,3) + a103*x*fastpow(z,3) + a203*fastpow(x,2)*fastpow(z,3) + a303*fastpow(x,3)*fastpow(z,3) + a013*y*fastpow(z,3) +
             a113*x*y*fastpow(z,3) + a213*fastpow(x,2)*y*fastpow(z,3) + a313*fastpow(x,3)*y*fastpow(z,3) + a023*fastpow(y,2)*fastpow(z,3) + a123*x*fastpow(y,2)*fastpow(z,3) +
             a223*fastpow(x,2)*fastpow(y,2)*fastpow(z,3) + a323*fastpow(x,3)*fastpow(y,2)*fastpow(z,3) + a033*fastpow(y,3)*fastpow(z,3) + a133*x*fastpow(y,3)*fastpow(z,3) +
             a233*fastpow(x,2)*fastpow(y,3)*fastpow(z,3) + a333*fastpow(x,3)*fastpow(y,3)*fastpow(z,3);
    }

  protected:
    void initCoeffsFromCellValues(const T cellValues[4][4][4]) {
      // This code has been generated from Mathematica solution to matching the values at the corners of the
      // domain, plus d/dx, d/dy, d/dz, d2/dxdy, d2/dxdz, d2/dydz, d3/dxdydz at each corner too.

      a000 = cellValues[1][1][1];
      a001 = (-cellValues[1][1][0] + cellValues[1][1][2]) / 2.;
      a002 = (2 * cellValues[1][1][0] - 5 * cellValues[1][1][1] + 4 * cellValues[1][1][2] - cellValues[1][1][3]) / 2.;
      a003 = (-cellValues[1][1][0] + 3 * cellValues[1][1][1] - 3 * cellValues[1][1][2] + cellValues[1][1][3]) / 2.;
      a010 = (-cellValues[1][0][1] + cellValues[1][2][1]) / 2.;
      a011 = (cellValues[1][0][0] - cellValues[1][0][2] - cellValues[1][2][0] + cellValues[1][2][2]) / 4.;
      a012 = (-2 * cellValues[1][0][0] + 5 * cellValues[1][0][1] - 4 * cellValues[1][0][2] + cellValues[1][0][3] +
              2 * cellValues[1][2][0] - 5 * cellValues[1][2][1] + 4 * cellValues[1][2][2] - cellValues[1][2][3]) / 4.;
      a013 = (cellValues[1][0][0] - 3 * cellValues[1][0][1] + 3 * cellValues[1][0][2] - cellValues[1][0][3] -
              cellValues[1][2][0] + 3 * cellValues[1][2][1] - 3 * cellValues[1][2][2] + cellValues[1][2][3]) / 4.;
      a020 = (2 * cellValues[1][0][1] - 5 * cellValues[1][1][1] + 4 * cellValues[1][2][1] - cellValues[1][3][1]) / 2.;
      a021 = (-2 * cellValues[1][0][0] + 2 * cellValues[1][0][2] + 5 * cellValues[1][1][0] - 5 * cellValues[1][1][2] -
              4 * cellValues[1][2][0] + 4 * cellValues[1][2][2] + cellValues[1][3][0] - cellValues[1][3][2]) / 4.;
      a022 = (4 * cellValues[1][0][0] - 10 * cellValues[1][0][1] + 8 * cellValues[1][0][2] - 2 * cellValues[1][0][3] -
              10 * cellValues[1][1][0] + 25 * cellValues[1][1][1] - 20 * cellValues[1][1][2] + 5 * cellValues[1][1][3] +
              8 * cellValues[1][2][0] - 20 * cellValues[1][2][1] + 16 * cellValues[1][2][2] - 4 * cellValues[1][2][3] -
              2 * cellValues[1][3][0] + 5 * cellValues[1][3][1] - 4 * cellValues[1][3][2] + cellValues[1][3][3]) / 4.;
      a023 = (-2 * cellValues[1][0][0] + 6 * cellValues[1][0][1] - 6 * cellValues[1][0][2] + 2 * cellValues[1][0][3] +
              5 * cellValues[1][1][0] - 15 * cellValues[1][1][1] + 15 * cellValues[1][1][2] - 5 * cellValues[1][1][3] -
              4 * cellValues[1][2][0] + 12 * cellValues[1][2][1] - 12 * cellValues[1][2][2] + 4 * cellValues[1][2][3] +
              cellValues[1][3][0] - 3 * cellValues[1][3][1] + 3 * cellValues[1][3][2] - cellValues[1][3][3]) / 4.;
      a030 = (-cellValues[1][0][1] + 3 * cellValues[1][1][1] - 3 * cellValues[1][2][1] + cellValues[1][3][1]) / 2.;
      a031 = (cellValues[1][0][0] - cellValues[1][0][2] - 3 * cellValues[1][1][0] + 3 * cellValues[1][1][2] +
              3 * cellValues[1][2][0] - 3 * cellValues[1][2][2] - cellValues[1][3][0] + cellValues[1][3][2]) / 4.;
      a032 = (-2 * cellValues[1][0][0] + 5 * cellValues[1][0][1] - 4 * cellValues[1][0][2] + cellValues[1][0][3] +
              6 * cellValues[1][1][0] - 15 * cellValues[1][1][1] + 12 * cellValues[1][1][2] - 3 * cellValues[1][1][3] -
              6 * cellValues[1][2][0] + 15 * cellValues[1][2][1] - 12 * cellValues[1][2][2] + 3 * cellValues[1][2][3] +
              2 * cellValues[1][3][0] - 5 * cellValues[1][3][1] + 4 * cellValues[1][3][2] - cellValues[1][3][3]) / 4.;
      a033 = (cellValues[1][0][0] - 3 * cellValues[1][0][1] + 3 * cellValues[1][0][2] - cellValues[1][0][3] -
              3 * cellValues[1][1][0] + 9 * cellValues[1][1][1] - 9 * cellValues[1][1][2] + 3 * cellValues[1][1][3] +
              3 * cellValues[1][2][0] - 9 * cellValues[1][2][1] + 9 * cellValues[1][2][2] - 3 * cellValues[1][2][3] -
              cellValues[1][3][0] + 3 * cellValues[1][3][1] - 3 * cellValues[1][3][2] + cellValues[1][3][3]) / 4.;
      a100 = (-cellValues[0][1][1] + cellValues[2][1][1]) / 2.;
      a101 = (cellValues[0][1][0] - cellValues[0][1][2] - cellValues[2][1][0] + cellValues[2][1][2]) / 4.;
      a102 = (-2 * cellValues[0][1][0] + 5 * cellValues[0][1][1] - 4 * cellValues[0][1][2] + cellValues[0][1][3] +
              2 * cellValues[2][1][0] - 5 * cellValues[2][1][1] + 4 * cellValues[2][1][2] - cellValues[2][1][3]) / 4.;
      a103 = (cellValues[0][1][0] - 3 * cellValues[0][1][1] + 3 * cellValues[0][1][2] - cellValues[0][1][3] -
              cellValues[2][1][0] + 3 * cellValues[2][1][1] - 3 * cellValues[2][1][2] + cellValues[2][1][3]) / 4.;
      a110 = (cellValues[0][0][1] - cellValues[0][2][1] - cellValues[2][0][1] + cellValues[2][2][1]) / 4.;
      a111 =
        (-cellValues[0][0][0] + cellValues[0][0][2] + cellValues[0][2][0] - cellValues[0][2][2] + cellValues[2][0][0] -
         cellValues[2][0][2] - cellValues[2][2][0] + cellValues[2][2][2]) / 8.;
      a112 = (2 * cellValues[0][0][0] - 5 * cellValues[0][0][1] + 4 * cellValues[0][0][2] - cellValues[0][0][3] -
              2 * cellValues[0][2][0] + 5 * cellValues[0][2][1] - 4 * cellValues[0][2][2] + cellValues[0][2][3] -
              2 * cellValues[2][0][0] + 5 * cellValues[2][0][1] - 4 * cellValues[2][0][2] + cellValues[2][0][3] +
              2 * cellValues[2][2][0] - 5 * cellValues[2][2][1] + 4 * cellValues[2][2][2] - cellValues[2][2][3]) / 8.;
      a113 = (-cellValues[0][0][0] + 3 * cellValues[0][0][1] - 3 * cellValues[0][0][2] + cellValues[0][0][3] +
              cellValues[0][2][0] - 3 * cellValues[0][2][1] + 3 * cellValues[0][2][2] - cellValues[0][2][3] +
              cellValues[2][0][0] - 3 * cellValues[2][0][1] + 3 * cellValues[2][0][2] - cellValues[2][0][3] -
              cellValues[2][2][0] + 3 * cellValues[2][2][1] - 3 * cellValues[2][2][2] + cellValues[2][2][3]) / 8.;
      a120 = (-2 * cellValues[0][0][1] + 5 * cellValues[0][1][1] - 4 * cellValues[0][2][1] + cellValues[0][3][1] +
              2 * cellValues[2][0][1] - 5 * cellValues[2][1][1] + 4 * cellValues[2][2][1] - cellValues[2][3][1]) / 4.;
      a121 = (2 * cellValues[0][0][0] - 2 * cellValues[0][0][2] - 5 * cellValues[0][1][0] + 5 * cellValues[0][1][2] +
              4 * cellValues[0][2][0] - 4 * cellValues[0][2][2] - cellValues[0][3][0] + cellValues[0][3][2] -
              2 * cellValues[2][0][0] + 2 * cellValues[2][0][2] + 5 * cellValues[2][1][0] - 5 * cellValues[2][1][2] -
              4 * cellValues[2][2][0] + 4 * cellValues[2][2][2] + cellValues[2][3][0] - cellValues[2][3][2]) / 8.;
      a122 = (-4 * cellValues[0][0][0] + 10 * cellValues[0][0][1] - 8 * cellValues[0][0][2] + 2 * cellValues[0][0][3] +
              10 * cellValues[0][1][0] - 25 * cellValues[0][1][1] + 20 * cellValues[0][1][2] - 5 * cellValues[0][1][3] -
              8 * cellValues[0][2][0] + 20 * cellValues[0][2][1] - 16 * cellValues[0][2][2] + 4 * cellValues[0][2][3] +
              2 * cellValues[0][3][0] - 5 * cellValues[0][3][1] + 4 * cellValues[0][3][2] - cellValues[0][3][3] +
              4 * cellValues[2][0][0] - 10 * cellValues[2][0][1] + 8 * cellValues[2][0][2] - 2 * cellValues[2][0][3] -
              10 * cellValues[2][1][0] + 25 * cellValues[2][1][1] - 20 * cellValues[2][1][2] + 5 * cellValues[2][1][3] +
              8 * cellValues[2][2][0] - 20 * cellValues[2][2][1] + 16 * cellValues[2][2][2] - 4 * cellValues[2][2][3] -
              2 * cellValues[2][3][0] + 5 * cellValues[2][3][1] - 4 * cellValues[2][3][2] + cellValues[2][3][3]) / 8.;
      a123 = (2 * cellValues[0][0][0] - 6 * cellValues[0][0][1] + 6 * cellValues[0][0][2] - 2 * cellValues[0][0][3] -
              5 * cellValues[0][1][0] + 15 * cellValues[0][1][1] - 15 * cellValues[0][1][2] + 5 * cellValues[0][1][3] +
              4 * cellValues[0][2][0] - 12 * cellValues[0][2][1] + 12 * cellValues[0][2][2] - 4 * cellValues[0][2][3] -
              cellValues[0][3][0] + 3 * cellValues[0][3][1] - 3 * cellValues[0][3][2] + cellValues[0][3][3] -
              2 * cellValues[2][0][0] + 6 * cellValues[2][0][1] - 6 * cellValues[2][0][2] + 2 * cellValues[2][0][3] +
              5 * cellValues[2][1][0] - 15 * cellValues[2][1][1] + 15 * cellValues[2][1][2] - 5 * cellValues[2][1][3] -
              4 * cellValues[2][2][0] + 12 * cellValues[2][2][1] - 12 * cellValues[2][2][2] + 4 * cellValues[2][2][3] +
              cellValues[2][3][0] - 3 * cellValues[2][3][1] + 3 * cellValues[2][3][2] - cellValues[2][3][3]) / 8.;
      a130 = (cellValues[0][0][1] - 3 * cellValues[0][1][1] + 3 * cellValues[0][2][1] - cellValues[0][3][1] -
              cellValues[2][0][1] + 3 * cellValues[2][1][1] - 3 * cellValues[2][2][1] + cellValues[2][3][1]) / 4.;
      a131 = (-cellValues[0][0][0] + cellValues[0][0][2] + 3 * cellValues[0][1][0] - 3 * cellValues[0][1][2] -
              3 * cellValues[0][2][0] + 3 * cellValues[0][2][2] + cellValues[0][3][0] - cellValues[0][3][2] +
              cellValues[2][0][0] - cellValues[2][0][2] - 3 * cellValues[2][1][0] + 3 * cellValues[2][1][2] +
              3 * cellValues[2][2][0] - 3 * cellValues[2][2][2] - cellValues[2][3][0] + cellValues[2][3][2]) / 8.;
      a132 = (2 * cellValues[0][0][0] - 5 * cellValues[0][0][1] + 4 * cellValues[0][0][2] - cellValues[0][0][3] -
              6 * cellValues[0][1][0] + 15 * cellValues[0][1][1] - 12 * cellValues[0][1][2] + 3 * cellValues[0][1][3] +
              6 * cellValues[0][2][0] - 15 * cellValues[0][2][1] + 12 * cellValues[0][2][2] - 3 * cellValues[0][2][3] -
              2 * cellValues[0][3][0] + 5 * cellValues[0][3][1] - 4 * cellValues[0][3][2] + cellValues[0][3][3] -
              2 * cellValues[2][0][0] + 5 * cellValues[2][0][1] - 4 * cellValues[2][0][2] + cellValues[2][0][3] +
              6 * cellValues[2][1][0] - 15 * cellValues[2][1][1] + 12 * cellValues[2][1][2] - 3 * cellValues[2][1][3] -
              6 * cellValues[2][2][0] + 15 * cellValues[2][2][1] - 12 * cellValues[2][2][2] + 3 * cellValues[2][2][3] +
              2 * cellValues[2][3][0] - 5 * cellValues[2][3][1] + 4 * cellValues[2][3][2] - cellValues[2][3][3]) / 8.;
      a133 = (-cellValues[0][0][0] + 3 * cellValues[0][0][1] - 3 * cellValues[0][0][2] + cellValues[0][0][3] +
              3 * cellValues[0][1][0] - 9 * cellValues[0][1][1] + 9 * cellValues[0][1][2] - 3 * cellValues[0][1][3] -
              3 * cellValues[0][2][0] + 9 * cellValues[0][2][1] - 9 * cellValues[0][2][2] + 3 * cellValues[0][2][3] +
              cellValues[0][3][0] - 3 * cellValues[0][3][1] + 3 * cellValues[0][3][2] - cellValues[0][3][3] +
              cellValues[2][0][0] - 3 * cellValues[2][0][1] + 3 * cellValues[2][0][2] - cellValues[2][0][3] -
              3 * cellValues[2][1][0] + 9 * cellValues[2][1][1] - 9 * cellValues[2][1][2] + 3 * cellValues[2][1][3] +
              3 * cellValues[2][2][0] - 9 * cellValues[2][2][1] + 9 * cellValues[2][2][2] - 3 * cellValues[2][2][3] -
              cellValues[2][3][0] + 3 * cellValues[2][3][1] - 3 * cellValues[2][3][2] + cellValues[2][3][3]) / 8.;
      a200 = (2 * cellValues[0][1][1] - 5 * cellValues[1][1][1] + 4 * cellValues[2][1][1] - cellValues[3][1][1]) / 2.;
      a201 = (-2 * cellValues[0][1][0] + 2 * cellValues[0][1][2] + 5 * cellValues[1][1][0] - 5 * cellValues[1][1][2] -
              4 * cellValues[2][1][0] + 4 * cellValues[2][1][2] + cellValues[3][1][0] - cellValues[3][1][2]) / 4.;
      a202 = (4 * cellValues[0][1][0] - 10 * cellValues[0][1][1] + 8 * cellValues[0][1][2] - 2 * cellValues[0][1][3] -
              10 * cellValues[1][1][0] + 25 * cellValues[1][1][1] - 20 * cellValues[1][1][2] + 5 * cellValues[1][1][3] +
              8 * cellValues[2][1][0] - 20 * cellValues[2][1][1] + 16 * cellValues[2][1][2] - 4 * cellValues[2][1][3] -
              2 * cellValues[3][1][0] + 5 * cellValues[3][1][1] - 4 * cellValues[3][1][2] + cellValues[3][1][3]) / 4.;
      a203 = (-2 * cellValues[0][1][0] + 6 * cellValues[0][1][1] - 6 * cellValues[0][1][2] + 2 * cellValues[0][1][3] +
              5 * cellValues[1][1][0] - 15 * cellValues[1][1][1] + 15 * cellValues[1][1][2] - 5 * cellValues[1][1][3] -
              4 * cellValues[2][1][0] + 12 * cellValues[2][1][1] - 12 * cellValues[2][1][2] + 4 * cellValues[2][1][3] +
              cellValues[3][1][0] - 3 * cellValues[3][1][1] + 3 * cellValues[3][1][2] - cellValues[3][1][3]) / 4.;
      a210 = (-2 * cellValues[0][0][1] + 2 * cellValues[0][2][1] + 5 * cellValues[1][0][1] - 5 * cellValues[1][2][1] -
              4 * cellValues[2][0][1] + 4 * cellValues[2][2][1] + cellValues[3][0][1] - cellValues[3][2][1]) / 4.;
      a211 = (2 * cellValues[0][0][0] - 2 * cellValues[0][0][2] - 2 * cellValues[0][2][0] + 2 * cellValues[0][2][2] -
              5 * cellValues[1][0][0] + 5 * cellValues[1][0][2] + 5 * cellValues[1][2][0] - 5 * cellValues[1][2][2] +
              4 * cellValues[2][0][0] - 4 * cellValues[2][0][2] - 4 * cellValues[2][2][0] + 4 * cellValues[2][2][2] -
              cellValues[3][0][0] + cellValues[3][0][2] + cellValues[3][2][0] - cellValues[3][2][2]) / 8.;
      a212 = (-4 * cellValues[0][0][0] + 10 * cellValues[0][0][1] - 8 * cellValues[0][0][2] + 2 * cellValues[0][0][3] +
              4 * cellValues[0][2][0] - 10 * cellValues[0][2][1] + 8 * cellValues[0][2][2] - 2 * cellValues[0][2][3] +
              10 * cellValues[1][0][0] - 25 * cellValues[1][0][1] + 20 * cellValues[1][0][2] - 5 * cellValues[1][0][3] -
              10 * cellValues[1][2][0] + 25 * cellValues[1][2][1] - 20 * cellValues[1][2][2] + 5 * cellValues[1][2][3] -
              8 * cellValues[2][0][0] + 20 * cellValues[2][0][1] - 16 * cellValues[2][0][2] + 4 * cellValues[2][0][3] +
              8 * cellValues[2][2][0] - 20 * cellValues[2][2][1] + 16 * cellValues[2][2][2] - 4 * cellValues[2][2][3] +
              2 * cellValues[3][0][0] - 5 * cellValues[3][0][1] + 4 * cellValues[3][0][2] - cellValues[3][0][3] -
              2 * cellValues[3][2][0] + 5 * cellValues[3][2][1] - 4 * cellValues[3][2][2] + cellValues[3][2][3]) / 8.;
      a213 = (2 * cellValues[0][0][0] - 6 * cellValues[0][0][1] + 6 * cellValues[0][0][2] - 2 * cellValues[0][0][3] -
              2 * cellValues[0][2][0] + 6 * cellValues[0][2][1] - 6 * cellValues[0][2][2] + 2 * cellValues[0][2][3] -
              5 * cellValues[1][0][0] + 15 * cellValues[1][0][1] - 15 * cellValues[1][0][2] + 5 * cellValues[1][0][3] +
              5 * cellValues[1][2][0] - 15 * cellValues[1][2][1] + 15 * cellValues[1][2][2] - 5 * cellValues[1][2][3] +
              4 * cellValues[2][0][0] - 12 * cellValues[2][0][1] + 12 * cellValues[2][0][2] - 4 * cellValues[2][0][3] -
              4 * cellValues[2][2][0] + 12 * cellValues[2][2][1] - 12 * cellValues[2][2][2] + 4 * cellValues[2][2][3] -
              cellValues[3][0][0] + 3 * cellValues[3][0][1] - 3 * cellValues[3][0][2] + cellValues[3][0][3] +
              cellValues[3][2][0] - 3 * cellValues[3][2][1] + 3 * cellValues[3][2][2] - cellValues[3][2][3]) / 8.;
      a220 = (4 * cellValues[0][0][1] - 10 * cellValues[0][1][1] + 8 * cellValues[0][2][1] - 2 * cellValues[0][3][1] -
              10 * cellValues[1][0][1] + 25 * cellValues[1][1][1] - 20 * cellValues[1][2][1] + 5 * cellValues[1][3][1] +
              8 * cellValues[2][0][1] - 20 * cellValues[2][1][1] + 16 * cellValues[2][2][1] - 4 * cellValues[2][3][1] -
              2 * cellValues[3][0][1] + 5 * cellValues[3][1][1] - 4 * cellValues[3][2][1] + cellValues[3][3][1]) / 4.;
      a221 = (-4 * cellValues[0][0][0] + 4 * cellValues[0][0][2] + 10 * cellValues[0][1][0] - 10 * cellValues[0][1][2] -
              8 * cellValues[0][2][0] + 8 * cellValues[0][2][2] + 2 * cellValues[0][3][0] - 2 * cellValues[0][3][2] +
              10 * cellValues[1][0][0] - 10 * cellValues[1][0][2] - 25 * cellValues[1][1][0] +
              25 * cellValues[1][1][2] + 20 * cellValues[1][2][0] - 20 * cellValues[1][2][2] - 5 * cellValues[1][3][0] +
              5 * cellValues[1][3][2] - 8 * cellValues[2][0][0] + 8 * cellValues[2][0][2] + 20 * cellValues[2][1][0] -
              20 * cellValues[2][1][2] - 16 * cellValues[2][2][0] + 16 * cellValues[2][2][2] + 4 * cellValues[2][3][0] -
              4 * cellValues[2][3][2] + 2 * cellValues[3][0][0] - 2 * cellValues[3][0][2] - 5 * cellValues[3][1][0] +
              5 * cellValues[3][1][2] + 4 * cellValues[3][2][0] - 4 * cellValues[3][2][2] - cellValues[3][3][0] +
              cellValues[3][3][2]) / 8.;
      a222 = (8 * cellValues[0][0][0] - 20 * cellValues[0][0][1] + 16 * cellValues[0][0][2] - 4 * cellValues[0][0][3] -
              20 * cellValues[0][1][0] + 50 * cellValues[0][1][1] - 40 * cellValues[0][1][2] +
              10 * cellValues[0][1][3] + 16 * cellValues[0][2][0] - 40 * cellValues[0][2][1] +
              32 * cellValues[0][2][2] - 8 * cellValues[0][2][3] - 4 * cellValues[0][3][0] + 10 * cellValues[0][3][1] -
              8 * cellValues[0][3][2] + 2 * cellValues[0][3][3] - 20 * cellValues[1][0][0] + 50 * cellValues[1][0][1] -
              40 * cellValues[1][0][2] + 10 * cellValues[1][0][3] + 50 * cellValues[1][1][0] -
              125 * cellValues[1][1][1] + 100 * cellValues[1][1][2] - 25 * cellValues[1][1][3] -
              40 * cellValues[1][2][0] + 100 * cellValues[1][2][1] - 80 * cellValues[1][2][2] +
              20 * cellValues[1][2][3] + 10 * cellValues[1][3][0] - 25 * cellValues[1][3][1] +
              20 * cellValues[1][3][2] - 5 * cellValues[1][3][3] + 16 * cellValues[2][0][0] - 40 * cellValues[2][0][1] +
              32 * cellValues[2][0][2] - 8 * cellValues[2][0][3] - 40 * cellValues[2][1][0] +
              100 * cellValues[2][1][1] - 80 * cellValues[2][1][2] + 20 * cellValues[2][1][3] +
              32 * cellValues[2][2][0] - 80 * cellValues[2][2][1] + 64 * cellValues[2][2][2] -
              16 * cellValues[2][2][3] - 8 * cellValues[2][3][0] + 20 * cellValues[2][3][1] - 16 * cellValues[2][3][2] +
              4 * cellValues[2][3][3] - 4 * cellValues[3][0][0] + 10 * cellValues[3][0][1] - 8 * cellValues[3][0][2] +
              2 * cellValues[3][0][3] + 10 * cellValues[3][1][0] - 25 * cellValues[3][1][1] + 20 * cellValues[3][1][2] -
              5 * cellValues[3][1][3] - 8 * cellValues[3][2][0] + 20 * cellValues[3][2][1] - 16 * cellValues[3][2][2] +
              4 * cellValues[3][2][3] + 2 * cellValues[3][3][0] - 5 * cellValues[3][3][1] + 4 * cellValues[3][3][2] -
              cellValues[3][3][3]) / 8.;
      a223 = (-4 * cellValues[0][0][0] + 12 * cellValues[0][0][1] - 12 * cellValues[0][0][2] + 4 * cellValues[0][0][3] +
              10 * cellValues[0][1][0] - 30 * cellValues[0][1][1] + 30 * cellValues[0][1][2] -
              10 * cellValues[0][1][3] - 8 * cellValues[0][2][0] + 24 * cellValues[0][2][1] - 24 * cellValues[0][2][2] +
              8 * cellValues[0][2][3] + 2 * cellValues[0][3][0] - 6 * cellValues[0][3][1] + 6 * cellValues[0][3][2] -
              2 * cellValues[0][3][3] + 10 * cellValues[1][0][0] - 30 * cellValues[1][0][1] + 30 * cellValues[1][0][2] -
              10 * cellValues[1][0][3] - 25 * cellValues[1][1][0] + 75 * cellValues[1][1][1] -
              75 * cellValues[1][1][2] + 25 * cellValues[1][1][3] + 20 * cellValues[1][2][0] -
              60 * cellValues[1][2][1] + 60 * cellValues[1][2][2] - 20 * cellValues[1][2][3] - 5 * cellValues[1][3][0] +
              15 * cellValues[1][3][1] - 15 * cellValues[1][3][2] + 5 * cellValues[1][3][3] - 8 * cellValues[2][0][0] +
              24 * cellValues[2][0][1] - 24 * cellValues[2][0][2] + 8 * cellValues[2][0][3] + 20 * cellValues[2][1][0] -
              60 * cellValues[2][1][1] + 60 * cellValues[2][1][2] - 20 * cellValues[2][1][3] -
              16 * cellValues[2][2][0] + 48 * cellValues[2][2][1] - 48 * cellValues[2][2][2] +
              16 * cellValues[2][2][3] + 4 * cellValues[2][3][0] - 12 * cellValues[2][3][1] + 12 * cellValues[2][3][2] -
              4 * cellValues[2][3][3] + 2 * cellValues[3][0][0] - 6 * cellValues[3][0][1] + 6 * cellValues[3][0][2] -
              2 * cellValues[3][0][3] - 5 * cellValues[3][1][0] + 15 * cellValues[3][1][1] - 15 * cellValues[3][1][2] +
              5 * cellValues[3][1][3] + 4 * cellValues[3][2][0] - 12 * cellValues[3][2][1] + 12 * cellValues[3][2][2] -
              4 * cellValues[3][2][3] - cellValues[3][3][0] + 3 * cellValues[3][3][1] - 3 * cellValues[3][3][2] +
              cellValues[3][3][3]) / 8.;
      a230 = (-2 * cellValues[0][0][1] + 6 * cellValues[0][1][1] - 6 * cellValues[0][2][1] + 2 * cellValues[0][3][1] +
              5 * cellValues[1][0][1] - 15 * cellValues[1][1][1] + 15 * cellValues[1][2][1] - 5 * cellValues[1][3][1] -
              4 * cellValues[2][0][1] + 12 * cellValues[2][1][1] - 12 * cellValues[2][2][1] + 4 * cellValues[2][3][1] +
              cellValues[3][0][1] - 3 * cellValues[3][1][1] + 3 * cellValues[3][2][1] - cellValues[3][3][1]) / 4.;
      a231 = (2 * cellValues[0][0][0] - 2 * cellValues[0][0][2] - 6 * cellValues[0][1][0] + 6 * cellValues[0][1][2] +
              6 * cellValues[0][2][0] - 6 * cellValues[0][2][2] - 2 * cellValues[0][3][0] + 2 * cellValues[0][3][2] -
              5 * cellValues[1][0][0] + 5 * cellValues[1][0][2] + 15 * cellValues[1][1][0] - 15 * cellValues[1][1][2] -
              15 * cellValues[1][2][0] + 15 * cellValues[1][2][2] + 5 * cellValues[1][3][0] - 5 * cellValues[1][3][2] +
              4 * cellValues[2][0][0] - 4 * cellValues[2][0][2] - 12 * cellValues[2][1][0] + 12 * cellValues[2][1][2] +
              12 * cellValues[2][2][0] - 12 * cellValues[2][2][2] - 4 * cellValues[2][3][0] + 4 * cellValues[2][3][2] -
              cellValues[3][0][0] + cellValues[3][0][2] + 3 * cellValues[3][1][0] - 3 * cellValues[3][1][2] -
              3 * cellValues[3][2][0] + 3 * cellValues[3][2][2] + cellValues[3][3][0] - cellValues[3][3][2]) / 8.;
      a232 = (-4 * cellValues[0][0][0] + 10 * cellValues[0][0][1] - 8 * cellValues[0][0][2] + 2 * cellValues[0][0][3] +
              12 * cellValues[0][1][0] - 30 * cellValues[0][1][1] + 24 * cellValues[0][1][2] - 6 * cellValues[0][1][3] -
              12 * cellValues[0][2][0] + 30 * cellValues[0][2][1] - 24 * cellValues[0][2][2] + 6 * cellValues[0][2][3] +
              4 * cellValues[0][3][0] - 10 * cellValues[0][3][1] + 8 * cellValues[0][3][2] - 2 * cellValues[0][3][3] +
              10 * cellValues[1][0][0] - 25 * cellValues[1][0][1] + 20 * cellValues[1][0][2] - 5 * cellValues[1][0][3] -
              30 * cellValues[1][1][0] + 75 * cellValues[1][1][1] - 60 * cellValues[1][1][2] +
              15 * cellValues[1][1][3] + 30 * cellValues[1][2][0] - 75 * cellValues[1][2][1] +
              60 * cellValues[1][2][2] - 15 * cellValues[1][2][3] - 10 * cellValues[1][3][0] +
              25 * cellValues[1][3][1] - 20 * cellValues[1][3][2] + 5 * cellValues[1][3][3] - 8 * cellValues[2][0][0] +
              20 * cellValues[2][0][1] - 16 * cellValues[2][0][2] + 4 * cellValues[2][0][3] + 24 * cellValues[2][1][0] -
              60 * cellValues[2][1][1] + 48 * cellValues[2][1][2] - 12 * cellValues[2][1][3] -
              24 * cellValues[2][2][0] + 60 * cellValues[2][2][1] - 48 * cellValues[2][2][2] +
              12 * cellValues[2][2][3] + 8 * cellValues[2][3][0] - 20 * cellValues[2][3][1] + 16 * cellValues[2][3][2] -
              4 * cellValues[2][3][3] + 2 * cellValues[3][0][0] - 5 * cellValues[3][0][1] + 4 * cellValues[3][0][2] -
              cellValues[3][0][3] - 6 * cellValues[3][1][0] + 15 * cellValues[3][1][1] - 12 * cellValues[3][1][2] +
              3 * cellValues[3][1][3] + 6 * cellValues[3][2][0] - 15 * cellValues[3][2][1] + 12 * cellValues[3][2][2] -
              3 * cellValues[3][2][3] - 2 * cellValues[3][3][0] + 5 * cellValues[3][3][1] - 4 * cellValues[3][3][2] +
              cellValues[3][3][3]) / 8.;
      a233 = (2 * cellValues[0][0][0] - 6 * cellValues[0][0][1] + 6 * cellValues[0][0][2] - 2 * cellValues[0][0][3] -
              6 * cellValues[0][1][0] + 18 * cellValues[0][1][1] - 18 * cellValues[0][1][2] + 6 * cellValues[0][1][3] +
              6 * cellValues[0][2][0] - 18 * cellValues[0][2][1] + 18 * cellValues[0][2][2] - 6 * cellValues[0][2][3] -
              2 * cellValues[0][3][0] + 6 * cellValues[0][3][1] - 6 * cellValues[0][3][2] + 2 * cellValues[0][3][3] -
              5 * cellValues[1][0][0] + 15 * cellValues[1][0][1] - 15 * cellValues[1][0][2] + 5 * cellValues[1][0][3] +
              15 * cellValues[1][1][0] - 45 * cellValues[1][1][1] + 45 * cellValues[1][1][2] -
              15 * cellValues[1][1][3] - 15 * cellValues[1][2][0] + 45 * cellValues[1][2][1] -
              45 * cellValues[1][2][2] + 15 * cellValues[1][2][3] + 5 * cellValues[1][3][0] - 15 * cellValues[1][3][1] +
              15 * cellValues[1][3][2] - 5 * cellValues[1][3][3] + 4 * cellValues[2][0][0] - 12 * cellValues[2][0][1] +
              12 * cellValues[2][0][2] - 4 * cellValues[2][0][3] - 12 * cellValues[2][1][0] + 36 * cellValues[2][1][1] -
              36 * cellValues[2][1][2] + 12 * cellValues[2][1][3] + 12 * cellValues[2][2][0] -
              36 * cellValues[2][2][1] + 36 * cellValues[2][2][2] - 12 * cellValues[2][2][3] - 4 * cellValues[2][3][0] +
              12 * cellValues[2][3][1] - 12 * cellValues[2][3][2] + 4 * cellValues[2][3][3] - cellValues[3][0][0] +
              3 * cellValues[3][0][1] - 3 * cellValues[3][0][2] + cellValues[3][0][3] + 3 * cellValues[3][1][0] -
              9 * cellValues[3][1][1] + 9 * cellValues[3][1][2] - 3 * cellValues[3][1][3] - 3 * cellValues[3][2][0] +
              9 * cellValues[3][2][1] - 9 * cellValues[3][2][2] + 3 * cellValues[3][2][3] + cellValues[3][3][0] -
              3 * cellValues[3][3][1] + 3 * cellValues[3][3][2] - cellValues[3][3][3]) / 8.;
      a300 = (-cellValues[0][1][1] + 3 * cellValues[1][1][1] - 3 * cellValues[2][1][1] + cellValues[3][1][1]) / 2.;
      a301 = (cellValues[0][1][0] - cellValues[0][1][2] - 3 * cellValues[1][1][0] + 3 * cellValues[1][1][2] +
              3 * cellValues[2][1][0] - 3 * cellValues[2][1][2] - cellValues[3][1][0] + cellValues[3][1][2]) / 4.;
      a302 = (-2 * cellValues[0][1][0] + 5 * cellValues[0][1][1] - 4 * cellValues[0][1][2] + cellValues[0][1][3] +
              6 * cellValues[1][1][0] - 15 * cellValues[1][1][1] + 12 * cellValues[1][1][2] - 3 * cellValues[1][1][3] -
              6 * cellValues[2][1][0] + 15 * cellValues[2][1][1] - 12 * cellValues[2][1][2] + 3 * cellValues[2][1][3] +
              2 * cellValues[3][1][0] - 5 * cellValues[3][1][1] + 4 * cellValues[3][1][2] - cellValues[3][1][3]) / 4.;
      a303 = (cellValues[0][1][0] - 3 * cellValues[0][1][1] + 3 * cellValues[0][1][2] - cellValues[0][1][3] -
              3 * cellValues[1][1][0] + 9 * cellValues[1][1][1] - 9 * cellValues[1][1][2] + 3 * cellValues[1][1][3] +
              3 * cellValues[2][1][0] - 9 * cellValues[2][1][1] + 9 * cellValues[2][1][2] - 3 * cellValues[2][1][3] -
              cellValues[3][1][0] + 3 * cellValues[3][1][1] - 3 * cellValues[3][1][2] + cellValues[3][1][3]) / 4.;
      a310 = (cellValues[0][0][1] - cellValues[0][2][1] - 3 * cellValues[1][0][1] + 3 * cellValues[1][2][1] +
              3 * cellValues[2][0][1] - 3 * cellValues[2][2][1] - cellValues[3][0][1] + cellValues[3][2][1]) / 4.;
      a311 = (-cellValues[0][0][0] + cellValues[0][0][2] + cellValues[0][2][0] - cellValues[0][2][2] +
              3 * cellValues[1][0][0] - 3 * cellValues[1][0][2] - 3 * cellValues[1][2][0] + 3 * cellValues[1][2][2] -
              3 * cellValues[2][0][0] + 3 * cellValues[2][0][2] + 3 * cellValues[2][2][0] - 3 * cellValues[2][2][2] +
              cellValues[3][0][0] - cellValues[3][0][2] - cellValues[3][2][0] + cellValues[3][2][2]) / 8.;
      a312 = (2 * cellValues[0][0][0] - 5 * cellValues[0][0][1] + 4 * cellValues[0][0][2] - cellValues[0][0][3] -
              2 * cellValues[0][2][0] + 5 * cellValues[0][2][1] - 4 * cellValues[0][2][2] + cellValues[0][2][3] -
              6 * cellValues[1][0][0] + 15 * cellValues[1][0][1] - 12 * cellValues[1][0][2] + 3 * cellValues[1][0][3] +
              6 * cellValues[1][2][0] - 15 * cellValues[1][2][1] + 12 * cellValues[1][2][2] - 3 * cellValues[1][2][3] +
              6 * cellValues[2][0][0] - 15 * cellValues[2][0][1] + 12 * cellValues[2][0][2] - 3 * cellValues[2][0][3] -
              6 * cellValues[2][2][0] + 15 * cellValues[2][2][1] - 12 * cellValues[2][2][2] + 3 * cellValues[2][2][3] -
              2 * cellValues[3][0][0] + 5 * cellValues[3][0][1] - 4 * cellValues[3][0][2] + cellValues[3][0][3] +
              2 * cellValues[3][2][0] - 5 * cellValues[3][2][1] + 4 * cellValues[3][2][2] - cellValues[3][2][3]) / 8.;
      a313 = (-cellValues[0][0][0] + 3 * cellValues[0][0][1] - 3 * cellValues[0][0][2] + cellValues[0][0][3] +
              cellValues[0][2][0] - 3 * cellValues[0][2][1] + 3 * cellValues[0][2][2] - cellValues[0][2][3] +
              3 * cellValues[1][0][0] - 9 * cellValues[1][0][1] + 9 * cellValues[1][0][2] - 3 * cellValues[1][0][3] -
              3 * cellValues[1][2][0] + 9 * cellValues[1][2][1] - 9 * cellValues[1][2][2] + 3 * cellValues[1][2][3] -
              3 * cellValues[2][0][0] + 9 * cellValues[2][0][1] - 9 * cellValues[2][0][2] + 3 * cellValues[2][0][3] +
              3 * cellValues[2][2][0] - 9 * cellValues[2][2][1] + 9 * cellValues[2][2][2] - 3 * cellValues[2][2][3] +
              cellValues[3][0][0] - 3 * cellValues[3][0][1] + 3 * cellValues[3][0][2] - cellValues[3][0][3] -
              cellValues[3][2][0] + 3 * cellValues[3][2][1] - 3 * cellValues[3][2][2] + cellValues[3][2][3]) / 8.;
      a320 = (-2 * cellValues[0][0][1] + 5 * cellValues[0][1][1] - 4 * cellValues[0][2][1] + cellValues[0][3][1] +
              6 * cellValues[1][0][1] - 15 * cellValues[1][1][1] + 12 * cellValues[1][2][1] - 3 * cellValues[1][3][1] -
              6 * cellValues[2][0][1] + 15 * cellValues[2][1][1] - 12 * cellValues[2][2][1] + 3 * cellValues[2][3][1] +
              2 * cellValues[3][0][1] - 5 * cellValues[3][1][1] + 4 * cellValues[3][2][1] - cellValues[3][3][1]) / 4.;
      a321 = (2 * cellValues[0][0][0] - 2 * cellValues[0][0][2] - 5 * cellValues[0][1][0] + 5 * cellValues[0][1][2] +
              4 * cellValues[0][2][0] - 4 * cellValues[0][2][2] - cellValues[0][3][0] + cellValues[0][3][2] -
              6 * cellValues[1][0][0] + 6 * cellValues[1][0][2] + 15 * cellValues[1][1][0] - 15 * cellValues[1][1][2] -
              12 * cellValues[1][2][0] + 12 * cellValues[1][2][2] + 3 * cellValues[1][3][0] - 3 * cellValues[1][3][2] +
              6 * cellValues[2][0][0] - 6 * cellValues[2][0][2] - 15 * cellValues[2][1][0] + 15 * cellValues[2][1][2] +
              12 * cellValues[2][2][0] - 12 * cellValues[2][2][2] - 3 * cellValues[2][3][0] + 3 * cellValues[2][3][2] -
              2 * cellValues[3][0][0] + 2 * cellValues[3][0][2] + 5 * cellValues[3][1][0] - 5 * cellValues[3][1][2] -
              4 * cellValues[3][2][0] + 4 * cellValues[3][2][2] + cellValues[3][3][0] - cellValues[3][3][2]) / 8.;
      a322 = (-4 * cellValues[0][0][0] + 10 * cellValues[0][0][1] - 8 * cellValues[0][0][2] + 2 * cellValues[0][0][3] +
              10 * cellValues[0][1][0] - 25 * cellValues[0][1][1] + 20 * cellValues[0][1][2] - 5 * cellValues[0][1][3] -
              8 * cellValues[0][2][0] + 20 * cellValues[0][2][1] - 16 * cellValues[0][2][2] + 4 * cellValues[0][2][3] +
              2 * cellValues[0][3][0] - 5 * cellValues[0][3][1] + 4 * cellValues[0][3][2] - cellValues[0][3][3] +
              12 * cellValues[1][0][0] - 30 * cellValues[1][0][1] + 24 * cellValues[1][0][2] - 6 * cellValues[1][0][3] -
              30 * cellValues[1][1][0] + 75 * cellValues[1][1][1] - 60 * cellValues[1][1][2] +
              15 * cellValues[1][1][3] + 24 * cellValues[1][2][0] - 60 * cellValues[1][2][1] +
              48 * cellValues[1][2][2] - 12 * cellValues[1][2][3] - 6 * cellValues[1][3][0] + 15 * cellValues[1][3][1] -
              12 * cellValues[1][3][2] + 3 * cellValues[1][3][3] - 12 * cellValues[2][0][0] + 30 * cellValues[2][0][1] -
              24 * cellValues[2][0][2] + 6 * cellValues[2][0][3] + 30 * cellValues[2][1][0] - 75 * cellValues[2][1][1] +
              60 * cellValues[2][1][2] - 15 * cellValues[2][1][3] - 24 * cellValues[2][2][0] +
              60 * cellValues[2][2][1] - 48 * cellValues[2][2][2] + 12 * cellValues[2][2][3] + 6 * cellValues[2][3][0] -
              15 * cellValues[2][3][1] + 12 * cellValues[2][3][2] - 3 * cellValues[2][3][3] + 4 * cellValues[3][0][0] -
              10 * cellValues[3][0][1] + 8 * cellValues[3][0][2] - 2 * cellValues[3][0][3] - 10 * cellValues[3][1][0] +
              25 * cellValues[3][1][1] - 20 * cellValues[3][1][2] + 5 * cellValues[3][1][3] + 8 * cellValues[3][2][0] -
              20 * cellValues[3][2][1] + 16 * cellValues[3][2][2] - 4 * cellValues[3][2][3] - 2 * cellValues[3][3][0] +
              5 * cellValues[3][3][1] - 4 * cellValues[3][3][2] + cellValues[3][3][3]) / 8.;
      a323 = (2 * cellValues[0][0][0] - 6 * cellValues[0][0][1] + 6 * cellValues[0][0][2] - 2 * cellValues[0][0][3] -
              5 * cellValues[0][1][0] + 15 * cellValues[0][1][1] - 15 * cellValues[0][1][2] + 5 * cellValues[0][1][3] +
              4 * cellValues[0][2][0] - 12 * cellValues[0][2][1] + 12 * cellValues[0][2][2] - 4 * cellValues[0][2][3] -
              cellValues[0][3][0] + 3 * cellValues[0][3][1] - 3 * cellValues[0][3][2] + cellValues[0][3][3] -
              6 * cellValues[1][0][0] + 18 * cellValues[1][0][1] - 18 * cellValues[1][0][2] + 6 * cellValues[1][0][3] +
              15 * cellValues[1][1][0] - 45 * cellValues[1][1][1] + 45 * cellValues[1][1][2] -
              15 * cellValues[1][1][3] - 12 * cellValues[1][2][0] + 36 * cellValues[1][2][1] -
              36 * cellValues[1][2][2] + 12 * cellValues[1][2][3] + 3 * cellValues[1][3][0] - 9 * cellValues[1][3][1] +
              9 * cellValues[1][3][2] - 3 * cellValues[1][3][3] + 6 * cellValues[2][0][0] - 18 * cellValues[2][0][1] +
              18 * cellValues[2][0][2] - 6 * cellValues[2][0][3] - 15 * cellValues[2][1][0] + 45 * cellValues[2][1][1] -
              45 * cellValues[2][1][2] + 15 * cellValues[2][1][3] + 12 * cellValues[2][2][0] -
              36 * cellValues[2][2][1] + 36 * cellValues[2][2][2] - 12 * cellValues[2][2][3] - 3 * cellValues[2][3][0] +
              9 * cellValues[2][3][1] - 9 * cellValues[2][3][2] + 3 * cellValues[2][3][3] - 2 * cellValues[3][0][0] +
              6 * cellValues[3][0][1] - 6 * cellValues[3][0][2] + 2 * cellValues[3][0][3] + 5 * cellValues[3][1][0] -
              15 * cellValues[3][1][1] + 15 * cellValues[3][1][2] - 5 * cellValues[3][1][3] - 4 * cellValues[3][2][0] +
              12 * cellValues[3][2][1] - 12 * cellValues[3][2][2] + 4 * cellValues[3][2][3] + cellValues[3][3][0] -
              3 * cellValues[3][3][1] + 3 * cellValues[3][3][2] - cellValues[3][3][3]) / 8.;
      a330 = (cellValues[0][0][1] - 3 * cellValues[0][1][1] + 3 * cellValues[0][2][1] - cellValues[0][3][1] -
              3 * cellValues[1][0][1] + 9 * cellValues[1][1][1] - 9 * cellValues[1][2][1] + 3 * cellValues[1][3][1] +
              3 * cellValues[2][0][1] - 9 * cellValues[2][1][1] + 9 * cellValues[2][2][1] - 3 * cellValues[2][3][1] -
              cellValues[3][0][1] + 3 * cellValues[3][1][1] - 3 * cellValues[3][2][1] + cellValues[3][3][1]) / 4.;
      a331 = (-cellValues[0][0][0] + cellValues[0][0][2] + 3 * cellValues[0][1][0] - 3 * cellValues[0][1][2] -
              3 * cellValues[0][2][0] + 3 * cellValues[0][2][2] + cellValues[0][3][0] - cellValues[0][3][2] +
              3 * cellValues[1][0][0] - 3 * cellValues[1][0][2] - 9 * cellValues[1][1][0] + 9 * cellValues[1][1][2] +
              9 * cellValues[1][2][0] - 9 * cellValues[1][2][2] - 3 * cellValues[1][3][0] + 3 * cellValues[1][3][2] -
              3 * cellValues[2][0][0] + 3 * cellValues[2][0][2] + 9 * cellValues[2][1][0] - 9 * cellValues[2][1][2] -
              9 * cellValues[2][2][0] + 9 * cellValues[2][2][2] + 3 * cellValues[2][3][0] - 3 * cellValues[2][3][2] +
              cellValues[3][0][0] - cellValues[3][0][2] - 3 * cellValues[3][1][0] + 3 * cellValues[3][1][2] +
              3 * cellValues[3][2][0] - 3 * cellValues[3][2][2] - cellValues[3][3][0] + cellValues[3][3][2]) / 8.;
      a332 = (2 * cellValues[0][0][0] - 5 * cellValues[0][0][1] + 4 * cellValues[0][0][2] - cellValues[0][0][3] -
              6 * cellValues[0][1][0] + 15 * cellValues[0][1][1] - 12 * cellValues[0][1][2] + 3 * cellValues[0][1][3] +
              6 * cellValues[0][2][0] - 15 * cellValues[0][2][1] + 12 * cellValues[0][2][2] - 3 * cellValues[0][2][3] -
              2 * cellValues[0][3][0] + 5 * cellValues[0][3][1] - 4 * cellValues[0][3][2] + cellValues[0][3][3] -
              6 * cellValues[1][0][0] + 15 * cellValues[1][0][1] - 12 * cellValues[1][0][2] + 3 * cellValues[1][0][3] +
              18 * cellValues[1][1][0] - 45 * cellValues[1][1][1] + 36 * cellValues[1][1][2] - 9 * cellValues[1][1][3] -
              18 * cellValues[1][2][0] + 45 * cellValues[1][2][1] - 36 * cellValues[1][2][2] + 9 * cellValues[1][2][3] +
              6 * cellValues[1][3][0] - 15 * cellValues[1][3][1] + 12 * cellValues[1][3][2] - 3 * cellValues[1][3][3] +
              6 * cellValues[2][0][0] - 15 * cellValues[2][0][1] + 12 * cellValues[2][0][2] - 3 * cellValues[2][0][3] -
              18 * cellValues[2][1][0] + 45 * cellValues[2][1][1] - 36 * cellValues[2][1][2] + 9 * cellValues[2][1][3] +
              18 * cellValues[2][2][0] - 45 * cellValues[2][2][1] + 36 * cellValues[2][2][2] - 9 * cellValues[2][2][3] -
              6 * cellValues[2][3][0] + 15 * cellValues[2][3][1] - 12 * cellValues[2][3][2] + 3 * cellValues[2][3][3] -
              2 * cellValues[3][0][0] + 5 * cellValues[3][0][1] - 4 * cellValues[3][0][2] + cellValues[3][0][3] +
              6 * cellValues[3][1][0] - 15 * cellValues[3][1][1] + 12 * cellValues[3][1][2] - 3 * cellValues[3][1][3] -
              6 * cellValues[3][2][0] + 15 * cellValues[3][2][1] - 12 * cellValues[3][2][2] + 3 * cellValues[3][2][3] +
              2 * cellValues[3][3][0] - 5 * cellValues[3][3][1] + 4 * cellValues[3][3][2] - cellValues[3][3][3]) / 8.;
      a333 = (-cellValues[0][0][0] + 3 * cellValues[0][0][1] - 3 * cellValues[0][0][2] + cellValues[0][0][3] +
              3 * cellValues[0][1][0] - 9 * cellValues[0][1][1] + 9 * cellValues[0][1][2] - 3 * cellValues[0][1][3] -
              3 * cellValues[0][2][0] + 9 * cellValues[0][2][1] - 9 * cellValues[0][2][2] + 3 * cellValues[0][2][3] +
              cellValues[0][3][0] - 3 * cellValues[0][3][1] + 3 * cellValues[0][3][2] - cellValues[0][3][3] +
              3 * cellValues[1][0][0] - 9 * cellValues[1][0][1] + 9 * cellValues[1][0][2] - 3 * cellValues[1][0][3] -
              9 * cellValues[1][1][0] + 27 * cellValues[1][1][1] - 27 * cellValues[1][1][2] + 9 * cellValues[1][1][3] +
              9 * cellValues[1][2][0] - 27 * cellValues[1][2][1] + 27 * cellValues[1][2][2] - 9 * cellValues[1][2][3] -
              3 * cellValues[1][3][0] + 9 * cellValues[1][3][1] - 9 * cellValues[1][3][2] + 3 * cellValues[1][3][3] -
              3 * cellValues[2][0][0] + 9 * cellValues[2][0][1] - 9 * cellValues[2][0][2] + 3 * cellValues[2][0][3] +
              9 * cellValues[2][1][0] - 27 * cellValues[2][1][1] + 27 * cellValues[2][1][2] - 9 * cellValues[2][1][3] -
              9 * cellValues[2][2][0] + 27 * cellValues[2][2][1] - 27 * cellValues[2][2][2] + 9 * cellValues[2][2][3] +
              3 * cellValues[2][3][0] - 9 * cellValues[2][3][1] + 9 * cellValues[2][3][2] - 3 * cellValues[2][3][3] +
              cellValues[3][0][0] - 3 * cellValues[3][0][1] + 3 * cellValues[3][0][2] - cellValues[3][0][3] -
              3 * cellValues[3][1][0] + 9 * cellValues[3][1][1] - 9 * cellValues[3][1][2] + 3 * cellValues[3][1][3] +
              3 * cellValues[3][2][0] - 9 * cellValues[3][2][1] + 9 * cellValues[3][2][2] - 3 * cellValues[3][2][3] -
              cellValues[3][3][0] + 3 * cellValues[3][3][1] - 3 * cellValues[3][3][2] + cellValues[3][3][3]) / 8.;
    }

  public:

    /* \brief Get the elements of the tranpose operation
     *
     * In the supporting notes, the matrix mapping a low-res onto an interpolated grid is \hat{P}^+.
     * This function helps construct \hat{P}^{+\dagger}, which is required during the process of constructing
     * multi-resolution covectors.
     *
     * Specifically, given a position in the unit cube, this routine fills in the weights associated with the
     * low-res pixels around the unit cube (those pixels that are actually used in constructing the interpolation
     * estimate). That means one can evaluate f(x) on every high-resolution pixel x, and then fval * f(x) are the
     * to be added to the downsampled field. Having done this for all points x in the high-resolution field, the
     * result is the tricubic-downsampled field, \hat{P}^{+\dagger} f.
     *
     * We can call this process 'conjugate deinterpolation', although note it is NOT an exact inverse since
     * \hat{P}^{+\dagger}\hat{P}^+ is not proportional to I.
     */
    static void getTransposeElementsForPosition(T x, T y, T z,  T fval[4][4][4])  {
      assert(x >= 0 && x < 1 && y >= 0 && y < 1 && z >= 0 && z < 1);
      fval[0][0][0] = -(fastpow(-1 + x, 2) * x * fastpow(-1 + y, 2) * y * fastpow(-1 + z, 2) * z) / 8.;
      fval[0][0][1] = (fastpow(-1 + x, 2) * x * fastpow(-1 + y, 2) * y * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[0][0][2] = -(fastpow(-1 + x, 2) * x * fastpow(-1 + y, 2) * y * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[0][0][3] = (fastpow(-1 + x, 2) * x * fastpow(-1 + y, 2) * y * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[0][1][0] = (fastpow(-1 + x, 2) * x * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[0][1][1] =
        -(fastpow(-1 + x, 2) * x * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[0][1][2] =
        (fastpow(-1 + x, 2) * x * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[0][1][3] = -(fastpow(-1 + x, 2) * x * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[0][2][0] = -(fastpow(-1 + x, 2) * x * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[0][2][1] =
        (fastpow(-1 + x, 2) * x * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[0][2][2] = -(fastpow(-1 + x, 2) * x * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[0][2][3] = (fastpow(-1 + x, 2) * x * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[0][3][0] = (fastpow(-1 + x, 2) * x * (-1 + y) * fastpow(y, 2) * fastpow(-1 + z, 2) * z) / 8.;
      fval[0][3][1] = -(fastpow(-1 + x, 2) * x * (-1 + y) * fastpow(y, 2) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[0][3][2] = (fastpow(-1 + x, 2) * x * (-1 + y) * fastpow(y, 2) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[0][3][3] = -(fastpow(-1 + x, 2) * x * (-1 + y) * fastpow(y, 2) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[1][0][0] = ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * fastpow(-1 + y, 2) * y * fastpow(-1 + z, 2) * z) / 8.;
      fval[1][0][1] =
        -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * fastpow(-1 + y, 2) * y * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[1][0][2] =
        ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * fastpow(-1 + y, 2) * y * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[1][0][3] = -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * fastpow(-1 + y, 2) * y * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[1][1][0] =
        -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[1][1][1] = ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) *
                       (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[1][1][2] = -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * z *
                        (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[1][1][3] =
        ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[1][2][0] =
        ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[1][2][1] = -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) *
                        (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[1][2][2] =
        ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) /
        8.;
      fval[1][2][3] =
        -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[1][3][0] = -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (-1 + y) * fastpow(y, 2) * fastpow(-1 + z, 2) * z) / 8.;
      fval[1][3][1] =
        ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (-1 + y) * fastpow(y, 2) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[1][3][2] =
        -((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (-1 + y) * fastpow(y, 2) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[1][3][3] = ((2 - 5 * fastpow(x, 2) + 3 * fastpow(x, 3)) * (-1 + y) * fastpow(y, 2) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[2][0][0] = -(x * (-1 - 4 * x + 3 * fastpow(x, 2)) * fastpow(-1 + y, 2) * y * fastpow(-1 + z, 2) * z) / 8.;
      fval[2][0][1] = (x * (-1 - 4 * x +
                            3 * fastpow(x, 2)) * fastpow(-1 + y, 2) * y * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[2][0][2] = -(x *
                        (-1 - 4 * x + 3 * fastpow(x, 2)) * fastpow(-1 + y, 2) * y * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[2][0][3] = (x * (-1 - 4 * x + 3 * fastpow(x, 2)) * fastpow(-1 + y, 2) * y * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[2][1][0] = (x * (-1 - 4 * x +
                            3 * fastpow(x, 2)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[2][1][1] = -(x * (-1 - 4 * x + 3 * fastpow(x, 2)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) *
                        (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[2][1][2] =
        (x * (-1 - 4 * x + 3 * fastpow(x, 2)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) /
        8.;
      fval[2][1][3] = -(x * (-1 - 4 * x +
                             3 * fastpow(x, 2)) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[2][2][0] = -(x *
                        (-1 - 4 * x + 3 * fastpow(x, 2)) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[2][2][1] = (x * (-1 - 4 * x + 3 * fastpow(x, 2)) * y *
                       (-1 - 4 * y + 3 * fastpow(y, 2)) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[2][2][2] = -(x * (-1 - 4 * x + 3 * fastpow(x, 2)) * y *
                        (-1 - 4 * y + 3 * fastpow(y, 2)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[2][2][3] = (x *
                       (-1 - 4 * x + 3 * fastpow(x, 2)) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[2][3][0] = (x * (-1 - 4 * x + 3 * fastpow(x, 2)) * (-1 + y) * fastpow(y, 2) * fastpow(-1 + z, 2) * z) / 8.;
      fval[2][3][1] = -(x * (-1 - 4 * x +
                             3 * fastpow(x, 2)) * (-1 + y) * fastpow(y, 2) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[2][3][2] = (x *
                       (-1 - 4 * x + 3 * fastpow(x, 2)) * (-1 + y) * fastpow(y, 2) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[2][3][3] = -(x * (-1 - 4 * x + 3 * fastpow(x, 2)) * (-1 + y) * fastpow(y, 2) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[3][0][0] = ((-1 + x) * fastpow(x, 2) * fastpow(-1 + y, 2) * y * fastpow(-1 + z, 2) * z) / 8.;
      fval[3][0][1] = -((-1 + x) * fastpow(x, 2) * fastpow(-1 + y, 2) * y * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[3][0][2] = ((-1 + x) * fastpow(x, 2) * fastpow(-1 + y, 2) * y * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[3][0][3] = -((-1 + x) * fastpow(x, 2) * fastpow(-1 + y, 2) * y * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[3][1][0] = -((-1 + x) * fastpow(x, 2) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[3][1][1] = ((-1 + x) *
                       fastpow(x, 2) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[3][1][2] = -((-1 + x) *
                        fastpow(x, 2) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[3][1][3] = ((-1 + x) * fastpow(x, 2) * (2 - 5 * fastpow(y, 2) + 3 * fastpow(y, 3)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[3][2][0] = ((-1 + x) * fastpow(x, 2) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * fastpow(-1 + z, 2) * z) / 8.;
      fval[3][2][1] = -((-1 + x) *
                        fastpow(x, 2) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[3][2][2] = ((-1 + x) *
                       fastpow(x, 2) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[3][2][3] = -((-1 + x) * fastpow(x, 2) * y * (-1 - 4 * y + 3 * fastpow(y, 2)) * (-1 + z) * fastpow(z, 2)) / 8.;
      fval[3][3][0] = -((-1 + x) * fastpow(x, 2) * (-1 + y) * fastpow(y, 2) * fastpow(-1 + z, 2) * z) / 8.;
      fval[3][3][1] = ((-1 + x) * fastpow(x, 2) * (-1 + y) * fastpow(y, 2) * (2 - 5 * fastpow(z, 2) + 3 * fastpow(z, 3))) / 8.;
      fval[3][3][2] = -((-1 + x) * fastpow(x, 2) * (-1 + y) * fastpow(y, 2) * z * (-1 - 4 * z + 3 * fastpow(z, 2))) / 8.;
      fval[3][3][3] = ((-1 + x) * fastpow(x, 2) * (-1 + y) * fastpow(y, 2) * (-1 + z) * fastpow(z, 2)) / 8.;
    }

  };
}
#endif //IC_TRICUBIC_HPP
