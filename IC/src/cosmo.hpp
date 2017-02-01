#ifndef _COSMO_HPP_INCLUDED
#define _COSMO_HPP_INCLUDED

#include <functional>
#include <climits>
#include "io.hpp"
#include "interpolation.hpp"
// #include "realspacecorrelation.hpp"

template<typename FloatType>
struct CosmologicalParameters {
  FloatType OmegaM0, OmegaLambda0, OmegaBaryons0, hubble, redshift;
  FloatType scalefactor, sigma8, ns, TCMB;
};

template<typename FloatType>
extern FloatType D(FloatType a, FloatType Om, FloatType Ol);

template<typename FloatType>
class CAMB {
protected:
  std::vector<FloatType> kcamb;
  std::vector<FloatType> Tcamb;
  Interpolator<FloatType> interpolator;
  FloatType amplitude;
  FloatType ns;
  mutable FloatType kcamb_max_in_file;

public:

  bool isUsable() const {
    return (kcamb.size() > 0);
  }

  void read(std::string incamb, const CosmologicalParameters<FloatType> &cosmology) {

    readLinesFromCambOutput(incamb);
    interpolator.initialise(kcamb, Tcamb);

    // a bit awkward that we have to copy this value:
    ns = cosmology.ns;

    calculateOverallNormalization(cosmology);


    /*
    realspace::RealSpaceGenerators<FloatType> obj(200.0,8192*2);

    cerr << "k=" << obj.generateKArray() << endl;
    cerr << "Pk=" << obj.generatePkArray(*this) << endl;
    */

  }

protected:
  void readLinesFromCambOutput(std::string incamb) {
    kcamb.clear();
    Tcamb.clear();

    const int c = 7; // number of columns in transfer function
    size_t j;

    std::vector<double> input;

    io::getBuffer(input, incamb);

    if (input.size() < c || input.size() % c != 0) {
      throw std::runtime_error("CAMB transfer file doesn't have a sensible number of rows and columns");
    }


    FloatType ap = input[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

    for (j = 0; j < input.size() / c; j++) {
      if (input[c * j] > 0) {
        // hard-coded to first two columns of CAMB file -
        kcamb.push_back(FloatType(input[c * j]));
        Tcamb.push_back(FloatType(input[c * j + 1]) / ap);
      }
      else continue;
    }

    size_t nCambLines = kcamb.size();

    // extend high-k range using Meszaros solution
    // This is a very naive approximation and a big warning will be issued if it's actually used

    FloatType Tcamb_f = Tcamb.back();
    kcamb_max_in_file = kcamb.back();
    FloatType keq = 0.01;
    while (kcamb.back() < 1000) {
      kcamb.push_back(kcamb.back() + 1.0);
      FloatType kratio = kcamb.back() / kcamb_max_in_file;
      Tcamb.push_back(Tcamb_f*pow(kratio,-2.0)*log(kcamb.back()/keq)/log(kcamb_max_in_file/keq));
    }



  }

public:

  FloatType operator()(FloatType k) const {
    FloatType linearTransfer;
    if (k != 0)
      linearTransfer = interpolator(k);
    else
      linearTransfer = 0.0;

    if(k>kcamb_max_in_file) {
      cerr << "WARNING: maximum k in CAMB input file is insufficient (" << kcamb_max_in_file << ")" << endl;
      cerr << "         extrapolating using naive Meszaros solution" << endl;
      kcamb_max_in_file = std::numeric_limits<FloatType>().max();
    }

    return amplitude * powf(k, ns) * linearTransfer * linearTransfer;
  }

  std::vector<FloatType> getPowerSpectrumForGrid(const Grid<FloatType> &grid) const {
    assert(kcamb.size() == Tcamb.size());

    FloatType norm = getPowerSpectrumNormalizationForGrid(grid);

    std::vector<FloatType> P(grid.size3);

    for (size_t i = 0; i < grid.size3; ++i) {
      FloatType k = grid.getFourierCellAbsK(i);
      P[i] = (*this)(k) * norm;
    }

    return P;

  }


protected:
  void calculateOverallNormalization(const CosmologicalParameters<FloatType> &cosmology) {
    FloatType growthFactor = D(cosmology.scalefactor, cosmology.OmegaM0, cosmology.OmegaLambda0);

    FloatType growthFactorNormalized = growthFactor / D(1., cosmology.OmegaM0, cosmology.OmegaLambda0);

    FloatType sigma8PreNormalization = sig(8., cosmology.ns);

    FloatType linearRenormFactor = (cosmology.sigma8 / sigma8PreNormalization) * growthFactorNormalized;

    amplitude = linearRenormFactor * linearRenormFactor;

  }

  FloatType getPowerSpectrumNormalizationForGrid(const Grid<FloatType> &grid) const {

    FloatType kw = 2. * M_PI / grid.boxsize;
    FloatType norm = kw * kw * kw / powf(2. * M_PI, 3.); //since kw=2pi/L, this is just 1/V_box

    return norm;
  }

public:


  FloatType sig(FloatType R, FloatType ns) {

    FloatType s = 0., k, t;

    FloatType amp = 9. / 2. / M_PI / M_PI;
    FloatType kmax = min(kcamb.back(), 200.0 / R);
    FloatType kmin = kcamb[0];

    FloatType dk = (kmax - kmin) / 50000.;
    for (k = kmin; k < kmax; k += dk) {

      t = interpolator(k);

      s += powf(k, ns + 2.) * ((sin(k * R) - k * R * cos(k * R)) / ((k * R) * (k * R) * (k * R))) *
           ((sin(k * R) - k * R * cos(k * R)) / ((k * R) * (k * R) * (k * R))) * t * t;

    }


    s = sqrt(s * amp * dk);
    return s;

  }

};


template<typename FloatType>
FloatType D(FloatType a, FloatType Om, FloatType Ol) {

  FloatType Hsq = Om / powf(a, 3.0) + (1. - Om - Ol) / a / a + Ol;
  FloatType d = 2.5 * a * Om / powf(a, 3.0) / Hsq / (powf(Om / Hsq / a / a / a, 4. / 7.) - Ol / Hsq +
                                                   (1. + 0.5 * Om / powf(a, 3.0) / Hsq) * (1. + 1. / 70. * Ol / Hsq));

  //simplify this...?

  return d;
}


template<typename FloatType>
void powsp(int n, std::complex<FloatType> *ft, const char *out, FloatType Boxlength) {

  int res = n;
  int nBins = 100;
  FloatType *inBin = new FloatType[nBins];
  FloatType *Gx = new FloatType[nBins];
  FloatType *kbin = new FloatType[nBins];
  FloatType kmax = M_PI / Boxlength * (FloatType) res, kmin = 2.0f * M_PI / (FloatType) Boxlength, dklog =
    log10(kmax / kmin) / nBins, kw = 2.0f * M_PI / (FloatType) Boxlength;

  int ix, iy, iz, idx, idx2;
  FloatType kfft;

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

        FloatType vabs = ft[idx].real() * ft[idx].real() + ft[idx].imag() * ft[idx].imag();

        int iix, iiy, iiz;

        if (ix > res / 2) iix = ix - res; else iix = ix;
        if (iy > res / 2) iiy = iy - res; else iiy = iy;
        if (iz > res / 2) iiz = iz - res; else iiz = iz;

        kfft = sqrt(iix * iix + iiy * iiy + iiz * iiz);
        FloatType k = kfft * kw;

        // correct for aliasing, formula from Jing (2005), ApJ 620, 559
        // assume isotropic aliasing (approx. true for k<kmax=knyquist)
        // this formula is for CIC interpolation scheme, which we use <- only needed for Powerspectrum

        FloatType JingCorr = (1.0f - 2.0f / 3.0f * sin(M_PI * k / kmax / 2.0f) * sin(M_PI * k / kmax / 2.0f));
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

  FloatType psnorm = powf(Boxlength / (2.0 * M_PI), 3.0);

  for (ix = 0; ix < nBins; ix++) {

    if (inBin[ix] > 0) {


      ofs << std::setw(16) << pow(10., log10(kmin) + dklog * (ix + 0.5))
      << std::setw(16) << kbin[ix] / inBin[ix]
      << std::setw(16) << (FloatType) (Gx[ix] / inBin[ix]) * psnorm
      << std::setw(16) << inBin[ix]
      << std::endl;

    }
  }
  ofs.close();

  free(inBin);
  free(kbin);
  free(Gx);

}

template<typename DataType, typename FloatType=strip_complex<DataType>>
void powsp_noJing(const Field<DataType> & field,
                  const std::vector<FloatType> &P0, const char *out, FloatType Boxlength) {

  int res = field.getGrid().size;
  int nBins = 100;
  std::vector<FloatType> inBin(nBins);
  std::vector<FloatType> kbin(nBins);
  std::vector<FloatType> Gx(nBins);
  std::vector<FloatType> Px(nBins);

  FloatType kmax = M_PI / Boxlength * (FloatType) res, kmin = 2.0f * M_PI / (FloatType) Boxlength, dklog =
    log10(kmax / kmin) / nBins, kw = 2.0f * M_PI / (FloatType) Boxlength;

  int ix, iy, iz, idx, idx2;
  FloatType kfft;

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

        auto fieldValue = field.evaluateFourierMode({ix,iy,iz});
        FloatType vabs = (fieldValue.real() * fieldValue.real() + fieldValue.imag() * fieldValue.imag());
        vabs*=fourier::getFourierCellWeight(field, idx);

        int iix, iiy, iiz;

        if (ix > res / 2) iix = ix - res; else iix = ix;
        if (iy > res / 2) iiy = iy - res; else iiy = iy;
        if (iz > res / 2) iiz = iz - res; else iiz = iz;

        kfft = sqrt(iix * iix + iiy * iiy + iiz * iiz);
        FloatType k = kfft * kw;

        //.. logarithmic spacing in k
        idx2 = (int) ((1.0f / dklog * log10(k / kmin)));

        if (k >= kmin && k < kmax) {

          Gx[idx2] += vabs / (FloatType) (res * res * res); //because FFT is now normalised with 1/sqrt(Ntot)
          Px[idx2] += P0[idx];
          kbin[idx2] += k;
          inBin[idx2]++;

        } else { continue; }

      }


  //... convert to physical units ...
  std::ofstream ofs(out);


  FloatType psnorm = powf(Boxlength / (2.0 * M_PI), 3.0);

  for (ix = 0; ix < nBins; ix++) {

    if (inBin[ix] > 0) {


      ofs << std::setw(16) << pow(10., log10(kmin) + dklog * (ix + 0.5))
      << std::setw(16) << kbin[ix] / inBin[ix]
      << std::setw(16) << (FloatType) (Px[ix] / inBin[ix]) * psnorm
      << std::setw(16) << (FloatType) (Gx[ix] / inBin[ix]) * psnorm
      << std::setw(16) << inBin[ix]
      << std::endl;

    }
  }
  ofs.close();


}

template<typename DataType, typename FloatType=strip_complex<DataType>>
DataType *poiss(DataType *out, DataType *in, int res,
                FloatType Boxlength, FloatType a, FloatType Om) {

  long i;
  FloatType prefac =
    3. / 2. * Om / a * 100. * 100. / (3. * 100000.) / (3. * 100000.); // =3/2 Om0/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)
  FloatType kw = 2.0f * M_PI / (FloatType) Boxlength, k_inv;

  int k1, k2, k3, kk1, kk2, kk3;

  for (k1 = 0; k1 < res; k1++) {
    for (k2 = 0; k2 < res; k2++) {
      for (k3 = 0; k3 < res; k3++) {
        // TODO: refactor this to use the provided routine for getting k coordinates in Grid class
        i = (k1 * res + k2) * (res) + k3;

        if (k1 > res / 2) kk1 = k1 - res; else kk1 = k1;
        if (k2 > res / 2) kk2 = k2 - res; else kk2 = k2;
        if (k3 > res / 2) kk3 = k3 - res; else kk3 = k3;

        k_inv = 1.0/((FloatType) (kk1 * kk1 + kk2 * kk2 + kk3 * kk3) * kw * kw);
        if(i==0)
          k_inv = 0;

        out[i] = -in[i] * prefac * k_inv;
      }
    }
  }

  return out;

}


/*
template<typename FloatType>
std::complex<FloatType> *rev_poiss(std::complex<FloatType> *out, std::complex<FloatType> *in, int res, FloatType Boxlength,
                                 FloatType a, FloatType Om) {

  long i;

  FloatType prefac =
    3. / 2. * Om / a * 100. * 100. / (3. * 100000.) / (3. * 100000.); // 3/2 Om/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)

  FloatType kw = 2.0f * M_PI / (FloatType) (Boxlength), k;

  int k1, k2, k3, kk1, kk2, kk3;

  for (k1 = 0; k1 < res; k1++) {
    for (k2 = 0; k2 < res; k2++) {
      for (k3 = 0; k3 < res; k3++) {
        i = (k1 * res + k2) * (res) + k3;

        if (k1 > res / 2) kk1 = k1 - res; else kk1 = k1;
        if (k2 > res / 2) kk2 = k2 - res; else kk2 = k2;
        if (k3 > res / 2) kk3 = k3 - res; else kk3 = k3;

        k = (FloatType) (kk1 * kk1 + kk2 * kk2 + kk3 * kk3) * kw * kw;

        out[i] = -k * in[i] / prefac;
      }
    }
  }

  out[0] = std::complex<FloatType>(0., 0.);

  return out;

}
 */

#endif
