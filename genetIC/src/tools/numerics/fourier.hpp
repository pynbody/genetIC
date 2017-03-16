#ifndef FFTH_INCLUDED
#define FFTH_INCLUDED

#ifdef FFTW3
// FFTW3 VERSION

#include <stdexcept>
#include <fftw3.h>

#ifdef _OPENMP

#include <omp.h>
#include "src/tools/util_functions.hpp"
#include "src/coordinate.hpp"

#endif

namespace grids {
  template<typename T>
  class Grid;
}

namespace fields {
  template<typename T, typename S>
  class Field;
}

namespace numerics {
  namespace fourier {
    void init_fftw_threads() {
#ifdef FFTW_THREADS
      if (fftw_init_threads() == 0)
        throw std::runtime_error("Cannot initialize FFTW threads");
#ifndef _OPENMP
      fftw_plan_with_nthreads(FFTW_THREADS);
#else
      fftw_plan_with_nthreads(omp_get_num_procs());
#endif
#endif
    }

    template<typename FloatType>
    void repackForReality(std::complex<FloatType> *ft_full, std::complex<FloatType> *ft_reduced, const int nx) {
      // for testing purposes only!
      const int nz = nx / 2 + 1;

      for (size_t x = 0; x < nx; x++)
        for (size_t y = 0; y < nx; y++)
          for (size_t z = 0; z < nz; z++)
            ft_reduced[z + nz * y + nz * nx * x] = ft_full[z + nx * y + nx * nx * x];

    }

    template<typename FloatType>
    void repackForReality(std::complex<FloatType> *ft, const int nx) {
      repackForReality(ft, ft, nx);
    }

    template<typename FloatType>
    void fft(std::complex<FloatType> *fto, std::complex<FloatType> *ftin,
             const unsigned int res, const int dir) {
      throw std::runtime_error("Sorry, the fourier transform has not been implemented for your specified precision");
      // you'll need to implement an alternative specialisation like the one below for the correct calls
      // see http://www.fftw.org/doc/Precision.html#Precision
    }

    template<typename FloatType>
    void fft(FloatType *fto, FloatType *ftin,
             const unsigned int res, const int dir) {
      throw std::runtime_error("Sorry, the fourier transform has not been implemented for your specified precision");
      // you'll need to implement an alternative specialisation like the one below for the correct calls
      // see http://www.fftw.org/doc/Precision.html#Precision
    }

    template<typename FloatType>
    std::vector<std::complex<FloatType>> fft_real_1d(std::vector<FloatType> in) {
      throw std::runtime_error("Sorry, the fourier transform has not been implemented for your specified precision");
      // you'll need to implement an alternative specialisation like the one below for the correct calls
      // see http://www.fftw.org/doc/Precision.html#Precision
    }

    template<typename FloatType>
    std::vector<FloatType> reverse_fft_real_1d(std::vector<std::complex<FloatType>> in) {
      throw std::runtime_error("Sorry, the fourier transform has not been implemented for your specified precision");
      // you'll need to implement an alternative specialisation like the one below for the correct calls
      // see http://www.fftw.org/doc/Precision.html#Precision
    }


    template<typename T>
    int getNyquistModeThatMustBeReal(const grids::Grid<T> &g);

    template<typename T>
    std::complex<T> getFourierCoefficient(const fields::Field<T, T> &field, int kx, int ky, int kz) {

      const grids::Grid<T> &g(field.getGrid());

      int abs_kx = std::abs(kx);
      int abs_ky = std::abs(ky);
      int abs_kz = std::abs(kz);

      int even_nyquist = getNyquistModeThatMustBeReal(g);

      size_t id_ppp = g.getCellIndex(Coordinate<int>(abs_kx, abs_ky, abs_kz));
      size_t id_ppm = g.getCellIndex(Coordinate<int>(abs_kx, abs_ky, -abs_kz));
      size_t id_pmp = g.getCellIndex(Coordinate<int>(abs_kx, -abs_ky, abs_kz));
      size_t id_pmm = g.getCellIndex(Coordinate<int>(abs_kx, -abs_ky, -abs_kz));
      size_t id_mpp = g.getCellIndex(Coordinate<int>(-abs_kx, abs_ky, abs_kz));
      size_t id_mpm = g.getCellIndex(Coordinate<int>(-abs_kx, abs_ky, -abs_kz));
      size_t id_mmp = g.getCellIndex(Coordinate<int>(-abs_kx, -abs_ky, abs_kz));
      size_t id_mmm = g.getCellIndex(Coordinate<int>(-abs_kx, -abs_ky, -abs_kz));

      T field_ppp = field[id_ppp];
      T field_ppm = field[id_ppm];
      T field_pmp = field[id_pmp];
      T field_pmm = field[id_pmm];
      T field_mpp = field[id_mpp];
      T field_mpm = field[id_mpm];
      T field_mmp = field[id_mmp];
      T field_mmm = field[id_mmm];

      // There is no imaginary part to the zero-frequency component

      if (abs_kx == 0 || abs_kx == even_nyquist) {
        field_mpp = 0;
        field_mpm = 0;
        field_mmp = 0;
        field_mmm = 0;
      }

      if (abs_ky == 0 || abs_ky == even_nyquist) {
        field_pmp = 0;
        field_pmm = 0;
        field_mmp = 0;
        field_mmm = 0;
      }

      if (abs_kz == 0 || abs_kz == even_nyquist) {
        field_ppm = 0;
        field_pmm = 0;
        field_mpm = 0;
        field_mmm = 0;
      }

      T real, imag;

      // Final linear combination into real and imag parts, derived in mathematica
      if (kx >= 0) {
        if (ky >= 0) {
          if (kz >= 0) {
            real = field_ppp - field_pmm - field_mpm - field_mmp;
            imag = field_ppm + field_pmp + field_mpp - field_mmm;
          } else {
            real = -field_mmp + field_mpm + field_pmm + field_ppp;
            imag = field_mmm + field_mpp + field_pmp - field_ppm;
          }
        } else {
          if (kz >= 0) {
            real = field_mmp - field_mpm + field_pmm + field_ppp;
            imag = field_mmm + field_mpp - field_pmp + field_ppm;
          } else {
            real = field_mmp + field_mpm - field_pmm + field_ppp;
            imag = -field_mmm + field_mpp - field_pmp - field_ppm;
          }
        }
      } else {
        if (ky >= 0) {
          if (kz >= 0) {
            real = field_mmp + field_mpm - field_pmm + field_ppp;
            imag = field_mmm - field_mpp + field_pmp + field_ppm;
          } else {
            real = field_mmp - field_mpm + field_pmm + field_ppp;
            imag = -field_mmm - field_mpp + field_pmp - field_ppm;
          }
        } else {
          if (kz >= 0) {
            real = -field_mmp + field_mpm + field_pmm + field_ppp;
            imag = -field_mmm - field_mpp - field_pmp + field_ppm;
          } else {
            real = -field_mmp - field_mpm - field_pmm + field_ppp;
            imag = field_mmm - field_mpp - field_pmp - field_ppm;
          }
        }
      }

      return std::complex<T>(real, imag);

    }

    template<typename T>
    int getNyquistModeThatMustBeReal(const grids::Grid<T> &g) {
      int even_nyquist;
      if (g.size % 2 == 0)
        even_nyquist = int(g.size) / 2;
      else
        even_nyquist = int(g.size) * 100; // arbitrary large number that will not be seen

      return even_nyquist;
    }


    template<typename T>
    void setFourierCoefficient(fields::Field<std::complex<T>, T> &field, int kx, int ky, int kz,
                               T realPart, T imagPart) {
      const grids::Grid<T> &g = field.getGrid();

      size_t id_k, id_negk;

      id_k = g.getCellIndex(Coordinate<int>(kx, ky, kz));
      id_negk = g.getCellIndex(Coordinate<int>(-kx, -ky, -kz));

      field[id_k] = std::complex<T>(realPart, imagPart);
      field[id_negk] = std::complex<T>(realPart, -imagPart);

    }


    template<typename T>
    void addFourierCoefficient(fields::Field<T, T> &field, int kx, int ky, int kz, T realPart, T imagPart) {
      if (kx < 0) {
        kx = -kx;
        ky = -ky;
        kz = -kz;
        imagPart = -imagPart;
      }

      const grids::Grid<T> &g(field.getGrid());

      int abs_kx = std::abs(kx);
      int abs_ky = std::abs(ky);
      int abs_kz = std::abs(kz);

      int even_nyquist = getNyquistModeThatMustBeReal(g);

      size_t id_ppp = g.getCellIndex(Coordinate<int>(abs_kx, abs_ky, abs_kz));
      size_t id_ppm = g.getCellIndex(Coordinate<int>(abs_kx, abs_ky, -abs_kz));
      size_t id_pmp = g.getCellIndex(Coordinate<int>(abs_kx, -abs_ky, abs_kz));
      size_t id_pmm = g.getCellIndex(Coordinate<int>(abs_kx, -abs_ky, -abs_kz));
      size_t id_mpp = g.getCellIndex(Coordinate<int>(-abs_kx, abs_ky, abs_kz));
      size_t id_mpm = g.getCellIndex(Coordinate<int>(-abs_kx, abs_ky, -abs_kz));
      size_t id_mmp = g.getCellIndex(Coordinate<int>(-abs_kx, -abs_ky, abs_kz));
      size_t id_mmm = g.getCellIndex(Coordinate<int>(-abs_kx, -abs_ky, -abs_kz));

      T ppp_weight, ppm_weight, pmp_weight, pmm_weight, mpp_weight, mpm_weight, mmp_weight, mmm_weight;
      ppp_weight = ppm_weight = pmp_weight = pmm_weight = mpp_weight = mpm_weight = mmp_weight = mmm_weight = 0.25;


      int adjustments = 0;

      if (abs_kx == 0 || abs_kx == even_nyquist) {
        mpp_weight = 0;
        mpm_weight = 0;
        mmp_weight = 0;
        mmm_weight = 0;
        ppp_weight *= 2;
        ppm_weight *= 2;
        pmp_weight *= 2;
        pmm_weight *= 2;
        adjustments++;
      }

      if (abs_ky == 0 || abs_ky == even_nyquist) {
        pmp_weight = 0;
        pmm_weight = 0;
        mmp_weight = 0;
        mmm_weight = 0;
        ppp_weight *= 2;
        ppm_weight *= 2;
        mpp_weight *= 2;
        mpm_weight *= 2;
        adjustments++;
      }

      if (abs_kz == 0 || abs_kz == even_nyquist) {
        ppm_weight = 0;
        pmm_weight = 0;
        mpm_weight = 0;
        mmm_weight = 0;
        if (adjustments < 2) {
          ppp_weight *= 2;
          pmp_weight *= 2;
          mpp_weight *= 2;
          mmp_weight *= 2;
        }
      }

      if (ky >= 0) {
        if (kz >= 0) {
          field[id_ppp] += ppp_weight * realPart;
          field[id_mpp] += imagPart * mpp_weight;
          field[id_ppm] += imagPart * ppm_weight;
          field[id_mpm] += -mpm_weight * realPart;
          field[id_pmp] += imagPart * pmp_weight;
          field[id_mmp] += -mmp_weight * realPart;
          field[id_pmm] += -pmm_weight * realPart;
          field[id_mmm] += -imagPart * mmm_weight;
        } else {
          field[id_ppp] += ppp_weight * realPart;
          field[id_ppm] += -imagPart * ppm_weight;
          field[id_pmp] += imagPart * pmp_weight;
          field[id_pmm] += pmm_weight * realPart;
          field[id_mpp] += imagPart * mpp_weight;
          field[id_mpm] += mpm_weight * realPart;
          field[id_mmp] += -mmp_weight * realPart;
          field[id_mmm] += imagPart * mmm_weight;
        }
      } else {
        if (kz >= 0) {
          field[id_ppp] += ppp_weight * realPart;
          field[id_ppm] += imagPart * ppm_weight;
          field[id_pmp] += -imagPart * pmp_weight;
          field[id_pmm] += pmm_weight * realPart;
          field[id_mpp] += imagPart * mpp_weight;
          field[id_mpm] += -mpm_weight * realPart;
          field[id_mmp] += mmp_weight * realPart;
          field[id_mmm] += imagPart * mmm_weight;
        } else {
          field[id_ppp] += ppp_weight * realPart;
          field[id_ppm] += -imagPart * ppm_weight;
          field[id_pmp] += -imagPart * pmp_weight;
          field[id_pmm] += -pmm_weight * realPart;
          field[id_mpp] += imagPart * mpp_weight;
          field[id_mpm] += mpm_weight * realPart;
          field[id_mmp] += mmp_weight * realPart;
          field[id_mmm] += -imagPart * mmm_weight;
        }
      }

    }

    template<typename T>
    std::complex<T> getFourierCoefficient(const fields::Field<std::complex<T>, T> &field, int kx, int ky, int kz) {
      return field[field.getGrid().getCellIndex(Coordinate<int>(kx, ky, kz))];
    }

    template<typename T>
    fields::Field<ensure_complex<T>, strip_complex<T>>
    getComplexFourierField(const fields::Field<T, strip_complex<T>> &field) {
      using complex_T = ensure_complex<T>;
      using underlying_T = strip_complex<T>;
      const grids::Grid<underlying_T> &g = field.getGrid();
      fields::Field<complex_T, underlying_T> out(g);

      std::cerr << "field-size=" << out.getDataVector().size() << std::endl;

      for (size_t i = 0; i < g.size3; ++i) {
        int kx, ky, kz;
        std::tie(kx, ky, kz) = g.getFourierCellCoordinate(i);
        out[i] = getFourierCoefficient(field, kx, ky, kz);
      }
      return out;
    };

    template<typename T>
    void setFourierCoefficient(fields::Field<T, T> &field, int kx, int ky, int kz,
                               strip_complex<T> realPart, strip_complex<T> imagPart) {
      // this may be *horribly* inefficient
      auto existingCoeff = getFourierCoefficient(field, kx, ky, kz);
      addFourierCoefficient(field, kx, ky, kz, realPart - existingCoeff.real(), imagPart - existingCoeff.imag());

      /*
      auto newCoeff = getFourierCoefficient(field, kx, ky, kz);
      assert (abs(realPart-newCoeff.real())<1e-10);
      assert (abs(imagPart-newCoeff.imag())<1e-10);
       */
    }


    template<typename T>
    void applyTransformationInFourierBasis(const fields::Field<std::complex<T>, T> &sourceField,
                                           std::function<std::tuple<std::complex<T>, std::complex<T>, std::complex<T>>(
                                             std::complex<T>, int, int, int)> transformation,
                                           fields::Field<std::complex<T>, T> &destField1,
                                           fields::Field<std::complex<T>, T> &destField2,
                                           fields::Field<std::complex<T>, T> &destField3) {

      const grids::Grid<T> &g = destField1.getGrid();
      assert (&g == &(sourceField.getGrid()));
      assert (&g == &(destField2.getGrid()));
      assert (&g == &(destField3.getGrid()));

      for (size_t i = 0; i < g.size3; ++i) {
        int kx, ky, kz;
        std::tie(kx, ky, kz) = g.getFourierCellCoordinate(i);
        std::tie(destField1[i], destField2[i], destField3[i]) = transformation(sourceField[i], kx, ky, kz);
      }
    }

    template<typename T>
    void applyTransformationInFourierBasis(const fields::Field<T, T> &sourceField,
                                           std::function<std::tuple<std::complex<T>, std::complex<T>, std::complex<T>>(
                                             std::complex<T>, int, int, int)> transformation,
                                           fields::Field<T, T> &destField1,
                                           fields::Field<T, T> &destField2,
                                           fields::Field<T, T> &destField3) {

      const grids::Grid<T> &g = destField1.getGrid();
      assert (&g == &(sourceField.getGrid()));
      assert (&g == &(destField2.getGrid()));
      assert (&g == &(destField3.getGrid()));

      std::fill(destField1.getDataVector().begin(), destField1.getDataVector().end(), 0);
      std::fill(destField2.getDataVector().begin(), destField2.getDataVector().end(), 0);
      std::fill(destField3.getDataVector().begin(), destField3.getDataVector().end(), 0);

      int size_by_2 = int(g.size) / 2;

      for (int kx = 0; kx <= size_by_2; kx++) {
        for (int ky = -size_by_2 + 1; ky <= size_by_2; ky++) {
          int lower_limit = -size_by_2 + 1;
          if ((kx == 0 || kx == size_by_2) && (ky == 0 || ky == size_by_2))
            lower_limit = 0;
          for (int kz = lower_limit; kz <= size_by_2; kz++) {
            std::complex<T> res1, res2, res3;
            auto sourceFieldValue = getFourierCoefficient(sourceField, kx, ky, kz);
            std::tie(res1, res2, res3) = transformation(sourceFieldValue, kx, ky, kz);
            setFourierCoefficient(destField1, kx, ky, kz, res1.real(), res1.imag());
            setFourierCoefficient(destField2, kx, ky, kz, res2.real(), res2.imag());
            setFourierCoefficient(destField3, kx, ky, kz, res3.real(), res3.imag());
          }
        }

      }
    }

    template<typename T>
    int getFourierCellWeight(const fields::Field<T, T> &field, size_t i) {
      const grids::Grid<T> &grid(field.getGrid());
      int kx, ky, kz;
      std::tie(kx, ky, kz) = grid.getFourierCellCoordinate(i);
      int even_nyquist = getNyquistModeThatMustBeReal(grid);

      int weight = 8;

      if (kx == 0 || kx == even_nyquist)
        weight >>= 1;
      if (ky == 0 || ky == even_nyquist)
        weight >>= 1;
      if (kz == 0 || kz == even_nyquist)
        weight >>= 1;

      return weight;
    }

    template<typename T>
    int getFourierCellWeight(const fields::Field<std::complex<T>, T> &field, size_t i) {
      return 1;
    }


    template<>
    void fft<double>(std::complex<double> *fto, std::complex<double> *ftin,
                     const unsigned int res, const int dir) {

      init_fftw_threads();

      fftw_plan plan;
      size_t i;
      double norm = pow(static_cast<double>(res), 1.5);
      size_t len = static_cast<size_t>(res * res);
      len *= res;

      if (dir == 1)
        plan = fftw_plan_dft_3d(res, res, res,
                                reinterpret_cast<fftw_complex *>(&ftin[0]),
                                reinterpret_cast<fftw_complex *>(&fto[0]),
                                FFTW_FORWARD, FFTW_ESTIMATE);

      else if (dir == -1)
        plan = fftw_plan_dft_3d(res, res, res,
                                reinterpret_cast<fftw_complex *>(&ftin[0]),
                                reinterpret_cast<fftw_complex *>(&fto[0]),
                                FFTW_BACKWARD, FFTW_ESTIMATE);


      else throw std::runtime_error("Incorrect direction parameter to fft");

      fftw_execute(plan);
      fftw_destroy_plan(plan);

#pragma omp parallel for schedule(static) private(i)
      for (i = 0; i < len; i++)
        fto[i] /= norm;


    }

    template<>
    void fft<double>(double *fto, double *ftin,
                     const unsigned int res, const int dir) {

      init_fftw_threads();

      fftw_plan plan;
      size_t i;
      double norm = pow(static_cast<double>(res), 1.5);
      size_t len = static_cast<size_t>(res * res);
      len *= res;

      if (dir == 1)
        plan = fftw_plan_r2r_3d(res, res, res,
                                reinterpret_cast<double *>(&ftin[0]),
                                reinterpret_cast<double *>(&fto[0]),
                                FFTW_R2HC, FFTW_R2HC, FFTW_R2HC,
                                FFTW_ESTIMATE);

      else if (dir == -1)
        plan = fftw_plan_r2r_3d(res, res, res,
                                reinterpret_cast<double *>(&ftin[0]),
                                reinterpret_cast<double *>(&fto[0]),
                                FFTW_HC2R, FFTW_HC2R, FFTW_HC2R,
                                FFTW_ESTIMATE);


      else throw std::runtime_error("Incorrect direction parameter to fft");

      fftw_execute(plan);
      fftw_destroy_plan(plan);

#pragma omp parallel for schedule(static) private(i)
      for (i = 0; i < len; i++)
        fto[i] /= norm;


    }

    std::vector<std::complex<double> > fft_1d(std::vector<std::complex<double>> input, const int dir) {

      init_fftw_threads();

      fftw_plan plan;


      std::vector<std::complex<double>> output;
      output.resize(input.size());
      plan = fftw_plan_dft_1d(input.size(),
                              reinterpret_cast<fftw_complex *>(&input[0]),
                              reinterpret_cast<fftw_complex *>(&output[0]),
                              dir, FFTW_ESTIMATE);


      fftw_execute(plan);
      fftw_destroy_plan(plan);

      output /= sqrt(double(input.size()));
      return output;
    }

    template<>
    std::vector<std::complex<double> > fft_real_1d<double>(std::vector<double> input) {

      init_fftw_threads();

      fftw_plan plan;


      std::vector<std::complex<double>> output;
      output.resize(input.size() / 2 + 1);
      plan = fftw_plan_dft_r2c_1d(input.size(),
                                  &input[0],
                                  reinterpret_cast<fftw_complex *>(&output[0]),
                                  FFTW_ESTIMATE);


      fftw_execute(plan);
      fftw_destroy_plan(plan);

      output /= sqrt(double(input.size()));
      return output;
    }

    template<>
    std::vector<double> reverse_fft_real_1d<double>(std::vector<std::complex<double>> input) {

      init_fftw_threads();

      fftw_plan plan;


      std::vector<double> output;
      if (input.back().imag() == 0)
        output.resize(2 * input.size() - 2);
      else
        output.resize(2 * input.size() - 1);

      plan = fftw_plan_dft_c2r_1d(output.size(),
                                  reinterpret_cast<fftw_complex *>(&input[0]),
                                  &output[0],
                                  FFTW_ESTIMATE);


      fftw_execute(plan);
      fftw_destroy_plan(plan);

      output /= sqrt(double(output.size()));
      return output;
    }

#else



    // FFTW2 VERSION

#ifndef DOUBLEPRECISION     /* default is single-precision */

#include <srfftw.h>

#ifdef HAVE_HDF5
            hid_t hdf_float = H5Tcopy (H5T_NATIVE_FLOAT);
            hid_t hdf_double = H5Tcopy (H5T_NATIVE_FLOAT);
#endif
        //#else
        //#if (DOUBLEPRECISION == 2)   /* mixed precision, do we want to implement this at all? */
        //typedef float   FloatType;
        //typedef double  MyDouble;
        //hid_t hdf_float = H5Tcopy (H5T_NATIVE_FLOAT);
        //hid_t hdf_double = H5Tcopy (H5T_NATIVE_DOUBLE);
#else                        /* everything double-precision */
#ifdef FFTW_TYPE_PREFIX
#include <drfftw.h>
#else
#include <rfftw.h>
#endif
#endif




    template<typename FloatType>
    void fft(std::complex<FloatType> *fto, std::complex<FloatType> *ftin, const int res, const int dir)
    { //works when fto and ftin and output are the same, but overwrites fto for the backwards transform so we have another temporary array there

        long i;

      if(dir==1){
        std::complex<FloatType> *fti = (std::complex<FloatType>*)calloc(res*res*res,sizeof(std::complex<FloatType>));
        for(i=0;i<res*res*res;i++){fti[i]=ftin[i]/sqrt((FloatType)(res*res*res));}
        //for(i=0;i<res*res*res;i++){fti[i]=ftin[i]/(FloatType)(res*res*res);}
        fftwnd_plan pf = fftw3d_create_plan( res, res, res, FFTW_FORWARD, FFTW_ESTIMATE );
        fftwnd_one(pf, reinterpret_cast<fftw_complex*>(&fti[0]), reinterpret_cast<fftw_complex*>(&fto[0]));
        fftwnd_destroy_plan(pf);
        free(fti);
        }

      else if(dir==-1){
        std::complex<FloatType> *fti = (std::complex<FloatType>*)calloc(res*res*res,sizeof(std::complex<FloatType>));
        //memcpy(fti,ftin,res*res*res*sizeof(fftw_complex));
        for(i=0;i<res*res*res;i++){fti[i]=ftin[i]/sqrt((FloatType)(res*res*res));}
        fftwnd_plan pb;
         pb = fftw3d_create_plan( res,res,res, FFTW_BACKWARD, FFTW_ESTIMATE );
        fftwnd_one(pb, reinterpret_cast<fftw_complex*>(&fti[0]), reinterpret_cast<fftw_complex*>(&fto[0]));
        fftwnd_destroy_plan(pb);
        free(fti);
          }

      else {std::cerr<<"wrong parameter for direction in fft"<<std::endl;}

    }

#endif

    unsigned int integerCubeRoot(unsigned long x) {
      int s;
      unsigned int y;
      unsigned long b;

      y = 0;
      for (s = 63; s >= 0; s -= 3) {
        y += y;
        b = 3 * y * ((unsigned long) y + 1) + 1;
        if ((x >> s) >= b) {
          x -= b << s;
          y++;
        }
      }
      return y;
    }

    template<typename T>
    void fft(std::vector<T> &fto, std::vector<T> &ftin, const int dir) {
      assert(fto.size() == ftin.size());
      size_t res = integerCubeRoot(fto.size());
      assert(res * res * res == fto.size());
      fft(fto.data(), ftin.data(), res, dir);
    }
  }
}

#endif // FFTH_INCLUDED
