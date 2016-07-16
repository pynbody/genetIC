#ifndef _COSMO_HPP_INCLUDED
#define _COSMO_HPP_INCLUDED

#include <functional>
#include "io.hpp"

template<typename MyFloat>
struct CosmologicalParameters {
  MyFloat OmegaM0, OmegaLambda0, OmegaBaryons0, hubble, redshift;
  MyFloat scalefactor, sigma8, ns, TCMB;
};

template<typename MyFloat>
extern MyFloat D(MyFloat a, MyFloat Om, MyFloat Ol);

template<typename MyFloat>
class CAMB {
protected:
  std::vector<MyFloat> kcamb;
  std::vector<MyFloat> Tcamb;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  void dealloc() {
    if (spline != NULL)
      gsl_spline_free(spline);
    if (acc != NULL)
      gsl_interp_accel_free(acc);
  }

public:

  CAMB() {
    spline = NULL;
    acc = NULL;
  }

  ~CAMB() {
    dealloc();
  }

  bool isUsable() {
    return (kcamb.size() > 0);
  }

  void read(std::string incamb) {
    kcamb.clear();
    Tcamb.clear();
    dealloc();
    const int c = 7; // number of columns in transfer function
    size_t j;

    std::vector<double> input;

    getBuffer(input, incamb);

    if (input.size() < c || input.size() % c != 0) {
      throw std::runtime_error("CAMB transfer file doesn't have a sensible number of rows and columns");
    }


    MyFloat ap = input[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

    for (j = 0; j < input.size() / c; j++) {
      if (input[c * j] > 0) {
        // hard-coded to first two columns of CAMB file -
        kcamb.push_back(MyFloat(input[c * j]));
        Tcamb.push_back(MyFloat(input[c * j + 1]) / ap);
      }
      else continue;
    }

    size_t nCambLines = kcamb.size();

    // extend high-k range using power law

    MyFloat gradient =
      log(Tcamb[nCambLines - 1] / Tcamb[nCambLines - 2]) / log(kcamb[nCambLines - 1] / kcamb[nCambLines - 2]);

    MyFloat Tcamb_f = Tcamb.back(), kcamb_f = kcamb.back();

    while (kcamb.back() < 1000) {
      kcamb.push_back(kcamb.back() + 1.0);
      Tcamb.push_back(exp(log(Tcamb_f) + gradient * (log(kcamb.back() / kcamb_f))));
    }


    assert(kcamb.size() == Tcamb.size());

    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, kcamb.size());
    gsl_spline_init(spline, kcamb.data(), Tcamb.data(), kcamb.size());

  }

  MyFloat getPowerSpectrumNormalizationForGrid(const CosmologicalParameters<MyFloat> &cosmology,
                                               const Grid<MyFloat> &grid) {
    MyFloat growthFactor = D(cosmology.scalefactor, cosmology.OmegaM0, cosmology.OmegaLambda0);

    MyFloat growthFactorNormalized = growthFactor / D(1., cosmology.OmegaM0, cosmology.OmegaLambda0);

    MyFloat sigma8PreNormalization = sig(8., cosmology.ns);

    MyFloat linearRenormFactor = (cosmology.sigma8 / sigma8PreNormalization) * growthFactorNormalized;

    MyFloat kw = 2. * M_PI / grid.boxsize;

    MyFloat amp = linearRenormFactor * linearRenormFactor;
    MyFloat norm = kw * kw * kw / powf(2. * M_PI, 3.); //since kw=2pi/L, this is just 1/V_box

    return norm * amp;
  }

  std::vector<MyFloat> getPowerSpectrumForGrid(const CosmologicalParameters<MyFloat> &cosmology,
                                               const Grid<MyFloat> &grid) {
    assert(kcamb.size() == Tcamb.size());

    MyFloat norm = getPowerSpectrumNormalizationForGrid(cosmology, grid);

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, kcamb.size());
    gsl_spline_init(spline, kcamb.data(), Tcamb.data(), kcamb.size());

    std::vector<MyFloat> P(grid.size3);

    for (size_t i = 0; i < grid.size3; ++i) {
      MyFloat k = grid.getFourierCellAbsK(i);
      MyFloat linearTransfer;
      if (k != 0)
        linearTransfer = gsl_spline_eval(spline, k, acc);
      else
        linearTransfer = 0.0;


      P[i] = powf(k, cosmology.ns) * linearTransfer * linearTransfer * norm;

    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return P;

  }


  MyFloat sig(MyFloat R, MyFloat ns) {

    MyFloat s = 0., k, t;

    MyFloat amp = 9. / 2. / M_PI / M_PI;
    MyFloat kmax = min(kcamb.back(), 200.0 / R);
    MyFloat kmin = kcamb[0];

    MyFloat dk = (kmax - kmin) / 50000.;
    for (k = kmin; k < kmax; k += dk) {

      t = gsl_spline_eval(spline, k, acc);

      s += powf(k, ns + 2.) * ((sin(k * R) - k * R * cos(k * R)) / ((k * R) * (k * R) * (k * R))) *
           ((sin(k * R) - k * R * cos(k * R)) / ((k * R) * (k * R) * (k * R))) * t * t;

    }


    s = sqrt(s * amp * dk);
    return s;

  }

  void applyTransfer(int res,
                     MyFloat kw, MyFloat ns, MyFloat norm_amp,
                     std::vector<std::complex<MyFloat>> &ft,
                     std::vector<std::complex<MyFloat>> &ftsc,
                     std::complex<MyFloat> *P,
                     Filter<MyFloat> &filter) {

    assert(kcamb.size() == Tcamb.size());

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, kcamb.size());
    gsl_spline_init(spline, kcamb.data(), Tcamb.data(), kcamb.size());

    MyFloat kk = 0.;
    long k1, k2, k3;
    long ii, jj, ll;
    long idk;

#pragma omp parallel for schedule(static) default(shared) private(k1,k2,k3,ii,jj,ll,idk,kk)
    for (k1 = -res / 2; k1 < res / 2; k1++) {
      if (k1 < 0) { ii = k1 + res; } else { ii = k1; }
      for (k2 = -res / 2; k2 < res / 2; k2++) {
        if (k2 < 0) { jj = k2 + res; } else { jj = k2; }
        for (k3 = -res / 2; k3 < res / 2; k3++) {
          if (k3 < 0) { ll = k3 + res; } else { ll = k3; }


          idk = (ii * (res) + jj) * (res) + ll;
          if (idk == 0) { continue; }
          kk = sqrt(k1 * k1 + k2 * k2 + k3 * k3) * kw;

          P[idk] = gsl_spline_eval(spline, kk, acc);
          P[idk] *= (P[idk] * std::complex<MyFloat>(powf(kk, ns)) * norm_amp);
          // P[idk]*=filter(kk);

          ftsc[idk] = sqrt(P[idk] * filter(kk)) * ft[idk];

        }
      }
    }


    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

  }
};


template<typename MyFloat>
MyFloat D(MyFloat a, MyFloat Om, MyFloat Ol) {

  MyFloat Hsq = Om / powf(a, 3.0) + (1. - Om - Ol) / a / a + Ol;
  MyFloat d = 2.5 * a * Om / powf(a, 3.0) / Hsq / (powf(Om / Hsq / a / a / a, 4. / 7.) - Ol / Hsq +
                                                   (1. + 0.5 * Om / powf(a, 3.0) / Hsq) * (1. + 1. / 70. * Ol / Hsq));

  //simplify this...?

  return d;
}


template<typename MyFloat>
void powsp(int n, std::complex<MyFloat> *ft, const char *out, MyFloat Boxlength) {

  int res = n;
  int nBins = 100;
  MyFloat *inBin = new MyFloat[nBins];
  MyFloat *Gx = new MyFloat[nBins];
  MyFloat *kbin = new MyFloat[nBins];
  MyFloat kmax = M_PI / Boxlength * (MyFloat) res, kmin = 2.0f * M_PI / (MyFloat) Boxlength, dklog =
    log10(kmax / kmin) / nBins, kw = 2.0f * M_PI / (MyFloat) Boxlength;

  int ix, iy, iz, idx, idx2;
  MyFloat kfft;

  for (ix = 0; ix < nBins; ix++) {
    inBin[ix] = 0;
    Gx[ix] = 0.0f;
    kbin[ix] = 0.0f;
  }


  for (ix = 0; ix < res; ix++)
    for (iy = 0; iy < res; iy++)
      for (iz = 0; iz < res; iz++) {
        idx = (ix * res + iy) * (res) + iz;

        // determine mode modulus

        MyFloat vabs = ft[idx].real() * ft[idx].real() + ft[idx].imag() * ft[idx].imag();

        int iix, iiy, iiz;

        if (ix > res / 2) iix = ix - res; else iix = ix;
        if (iy > res / 2) iiy = iy - res; else iiy = iy;
        if (iz > res / 2) iiz = iz - res; else iiz = iz;

        kfft = sqrt(iix * iix + iiy * iiy + iiz * iiz);
        MyFloat k = kfft * kw;

        // correct for aliasing, formula from Jing (2005), ApJ 620, 559
        // assume isotropic aliasing (approx. true for k<kmax=knyquist)
        // this formula is for CIC interpolation scheme, which we use <- only needed for Powerspectrum

        MyFloat JingCorr = (1.0f - 2.0f / 3.0f * sin(M_PI * k / kmax / 2.0f) * sin(M_PI * k / kmax / 2.0f));
        vabs /= JingCorr;

        //.. logarithmic spacing in k
        idx2 = (int) ((1.0f / dklog * log10(k / kmin))); //cout<<idx2<<endl;

        if (k >= kmin && k < kmax) {


          Gx[idx2] += vabs;
          kbin[idx2] += k;
          inBin[idx2]++;

        } else { continue; }

      }


  //... convert to physical units ...
  std::ofstream ofs(out);

  //definition of powerspectrum brings (2pi)^-3, FT+conversion to physical units brings sqrt(Box^3/N^6) per delta1, where ps22~d2*d2~d1*d1*d1*d1 -> (Box^3/N^6)^2

  MyFloat psnorm = powf(Boxlength / (2.0 * M_PI), 3.0);

  for (ix = 0; ix < nBins; ix++) {

    if (inBin[ix] > 0) {


      ofs << std::setw(16) << pow(10., log10(kmin) + dklog * (ix + 0.5))
      << std::setw(16) << kbin[ix] / inBin[ix]
      << std::setw(16) << (MyFloat) (Gx[ix] / inBin[ix]) * psnorm
      << std::setw(16) << inBin[ix]
      << std::endl;

    }
  }
  ofs.close();

  free(inBin);
  free(kbin);
  free(Gx);

}

template<typename MyFloat>
void powsp_noJing(int n, const std::vector<std::complex<MyFloat>> &ft,
                  const std::vector<MyFloat> &P0, const char *out, MyFloat Boxlength) {

  int res = n;
  int nBins = 100;
  std::vector<MyFloat> inBin(nBins);
  std::vector<MyFloat> kbin(nBins);
  std::vector<MyFloat> Gx(nBins);
  std::vector<MyFloat> Px(nBins);

  MyFloat kmax = M_PI / Boxlength * (MyFloat) res, kmin = 2.0f * M_PI / (MyFloat) Boxlength, dklog =
    log10(kmax / kmin) / nBins, kw = 2.0f * M_PI / (MyFloat) Boxlength;

  int ix, iy, iz, idx, idx2;
  MyFloat kfft;

  for (ix = 0; ix < nBins; ix++) {
    inBin[ix] = 0;
    Gx[ix] = 0.0f;
    kbin[ix] = 0.0f;
  }


  for (ix = 0; ix < res; ix++)
    for (iy = 0; iy < res; iy++)
      for (iz = 0; iz < res; iz++) {
        idx = (ix * res + iy) * (res) + iz;

        // determine mode modulus

        MyFloat vabs = ft[idx].real() * ft[idx].real() + ft[idx].imag() * ft[idx].imag();

        int iix, iiy, iiz;

        if (ix > res / 2) iix = ix - res; else iix = ix;
        if (iy > res / 2) iiy = iy - res; else iiy = iy;
        if (iz > res / 2) iiz = iz - res; else iiz = iz;

        kfft = sqrt(iix * iix + iiy * iiy + iiz * iiz);
        MyFloat k = kfft * kw;

        //.. logarithmic spacing in k
        idx2 = (int) ((1.0f / dklog * log10(k / kmin)));

        if (k >= kmin && k < kmax) {

          Gx[idx2] += vabs / (MyFloat) (res * res * res); //because FFT is now normalised with 1/sqrt(Ntot)
          Px[idx2] += P0[idx];
          kbin[idx2] += k;
          inBin[idx2]++;

        } else { continue; }

      }


  //... convert to physical units ...
  std::ofstream ofs(out);


  MyFloat psnorm = powf(Boxlength / (2.0 * M_PI), 3.0);

  for (ix = 0; ix < nBins; ix++) {

    if (inBin[ix] > 0) {


      ofs << std::setw(16) << pow(10., log10(kmin) + dklog * (ix + 0.5))
      << std::setw(16) << kbin[ix] / inBin[ix]
      << std::setw(16) << (MyFloat) (Px[ix] / inBin[ix]) * psnorm
      << std::setw(16) << (MyFloat) (Gx[ix] / inBin[ix]) * psnorm
      << std::setw(16) << inBin[ix]
      << std::endl;

    }
  }
  ofs.close();


}

template<typename MyFloat>
std::complex<MyFloat> *poiss(std::complex<MyFloat> *out, std::complex<MyFloat> *in, int res, MyFloat Boxlength,
                             MyFloat a, MyFloat Om) {

  long i;
  MyFloat prefac =
    3. / 2. * Om / a * 100. * 100. / (3. * 100000.) / (3. * 100000.); // =3/2 Om0/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)
  MyFloat kw = 2.0f * M_PI / (MyFloat) Boxlength, k;

  int k1, k2, k3, kk1, kk2, kk3;

  for (k1 = 0; k1 < res; k1++) {
    for (k2 = 0; k2 < res; k2++) {
      for (k3 = 0; k3 < res; k3++) {
        i = (k1 * res + k2) * (res) + k3;

        if (k1 > res / 2) kk1 = k1 - res; else kk1 = k1;
        if (k2 > res / 2) kk2 = k2 - res; else kk2 = k2;
        if (k3 > res / 2) kk3 = k3 - res; else kk3 = k3;

        k = (MyFloat) (kk1 * kk1 + kk2 * kk2 + kk3 * kk3) * kw * kw;

        out[i] = -in[i] * prefac / k;
      }
    }
  }

  out[0] = std::complex<MyFloat>(0., 0.);

  return out;

}

template<typename MyFloat>
std::complex<MyFloat> *rev_poiss(std::complex<MyFloat> *out, std::complex<MyFloat> *in, int res, MyFloat Boxlength,
                                 MyFloat a, MyFloat Om) {

  long i;

  MyFloat prefac =
    3. / 2. * Om / a * 100. * 100. / (3. * 100000.) / (3. * 100000.); // 3/2 Om/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)

  MyFloat kw = 2.0f * M_PI / (MyFloat) (Boxlength), k;

  int k1, k2, k3, kk1, kk2, kk3;

  for (k1 = 0; k1 < res; k1++) {
    for (k2 = 0; k2 < res; k2++) {
      for (k3 = 0; k3 < res; k3++) {
        i = (k1 * res + k2) * (res) + k3;

        if (k1 > res / 2) kk1 = k1 - res; else kk1 = k1;
        if (k2 > res / 2) kk2 = k2 - res; else kk2 = k2;
        if (k3 > res / 2) kk3 = k3 - res; else kk3 = k3;

        k = (MyFloat) (kk1 * kk1 + kk2 * kk2 + kk3 * kk3) * kw * kw;

        out[i] = -k * in[i] / prefac;
      }
    }
  }

  out[0] = std::complex<MyFloat>(0., 0.);

  return out;

}

#endif
