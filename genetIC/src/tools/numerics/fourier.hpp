#ifndef FFTH_INCLUDED
#define FFTH_INCLUDED


#include <stdexcept>
#include <fftw3.h>

#ifdef _OPENMP

#include <omp.h>

#endif

#include <iostream>
#include <src/simulation/field/field.hpp>

#include "src/simulation/coordinate.hpp"
#include "src/tools/data_types/complex.hpp"
#include "src/tools/numerics/vectormath.hpp"
#include "src/simulation/grid/grid.hpp"
#include "src/simulation/field/field.hpp"

namespace tools {
  namespace numerics {
    namespace fourier {
      bool fftwThreadsInitialised = false;

      void initialise() {
        if (fftwThreadsInitialised)
          return;

#ifdef FFTW_THREADS
        if (fftw_init_threads() == 0)
          throw std::runtime_error("Cannot initialize FFTW threads");
#ifndef _OPENMP
        fftw_plan_with_nthreads(FFTW_THREADS);
	std::cerr << "Note: " << FFTW_THREADS << " FFTW Threads were initialised" << std::endl;
#else
        fftw_plan_with_nthreads(omp_get_max_threads());
        std::cerr << "Note: " << omp_get_max_threads() << " FFTW Threads (determined by OpenMP) were initialised"
                  << std::endl;
#endif
#else
	std::cerr << "Note: FFTW Threads are not enabled" << std::endl;
#endif
        fftwThreadsInitialised = true;
      }

      template<>
      class FieldFourierManager<double> {
      protected:
        using T=double;
        fields::Field<double, double> &field;
        const grids::Grid<double> &grid;
        int size;
        size_t compressed_size;

        auto getRealCoeffLocationAndConjugation(int kx, int ky, int kz) {
          bool conjugate = false;
          if(kz<0) {
            conjugate = true;
            kx = -kx;
            ky = -ky;
            kz = -kz;
          }
          if(ky<0) {
            ky = size+ky;
          }
          if(kx<0) {
            kx = size+kx;
          }
          size_t logical_index = kz + compressed_size*ky + compressed_size*size_t(size*kx);
          size_t index_re = 2 * logical_index;
          assert(index_re+1<field.getDataVector().size());

          return std::make_tuple(conjugate,index_re);
        }

        void mirrorFourierModesForFFTWRealTransform() {
          // for conceptual ease, the FFTW real transformations contain some duplicate data
          // see  http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html
          //
          // we always address positive fourier modes where there is an ambiguity, but it is
          // not clear that FFTW does the same, so before a transform back to real space we need
          // to fill in the missing modes.
          size_t sourceLocation, destLocation;
          bool conj;

          // Copy missing data from kz=0, kx>0 to kz=0, kx<0
          // Conceptually: field [kx,ky,0] = conj(field[-kx,-ky,0])
          for(int ky=-size/2; ky<size/2+1; ++ky) {
            for(int kx=1; kx<size/2+1; ++kx) {
              std::tie(conj, sourceLocation) = getRealCoeffLocationAndConjugation(kx, ky, 0);
              std::tie(conj, destLocation) = getRealCoeffLocationAndConjugation(-kx, -ky, 0);
              field[destLocation] = field[sourceLocation];
              field[destLocation+1] = -field[sourceLocation+1];
            }
          }
          // Copy missing data from kz=0, ky>0, kx=0 to kz=0, ky<0, kx=0
          // Conceptually: field[0,ky,0] = conj(field[0,-ky,0])
          for(int ky=1; ky<size/2+1; ++ky) {
            std::tie(conj, sourceLocation) = getRealCoeffLocationAndConjugation(0, ky, 0);
            std::tie(conj, destLocation) = getRealCoeffLocationAndConjugation(0, -ky, 0);
            field[destLocation] = field[sourceLocation];
            field[destLocation+1] = -field[sourceLocation+1];
          }
        }

        void padForFFTWRealTransform() {

          // convert from our internal array format (which is just a standard column-major rep)
          // to the FFTW real-packed version which has padding every grid_size
          //
          // see http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html
          size_t source_range_end = grid.size3;
          size_t padding_amount = compressed_size*2-grid.size;
          size_t target_range_end = grid.size*grid.size*2*compressed_size-padding_amount;
          auto & data = field.getDataVector();

          while(source_range_end>grid.size) {
            size_t target_range_start = target_range_end - grid.size;
            size_t source_range_start = source_range_end - grid.size;
            std::copy_backward(&data[source_range_start], &data[source_range_end], &data[target_range_end]);
            target_range_end = target_range_start - padding_amount;
            source_range_end = source_range_start;
          }


        }

        void unpadAfterFFTWRealTransform() {
          // reverse of padForFFTWRealTransform

          size_t source_range_start = 0;
          size_t target_range_start = 0;
          size_t source_max = getRequiredDataSize();

          size_t padding_amount = compressed_size*2-grid.size;
          auto & data = field.getDataVector();

          while(source_range_start<source_max) {
            size_t source_range_end = source_range_start + grid.size;
            std::copy(&data[source_range_start], &data[source_range_end], &data[target_range_start]);
            source_range_start = source_range_end + padding_amount;
            target_range_start+=grid.size;
          }

        }
      public:
        FieldFourierManager(fields::Field<double, double> &field) : field(field), grid(field.getGrid()) {
          size = static_cast<int>(grid.size);
          compressed_size = grid.size/2+1;
        }

        void setFourierCoefficient(int kx, int ky, int kz, const std::complex<T> &val)  {
          bool conj;
          size_t index_re;

          std::tie(conj, index_re) = getRealCoeffLocationAndConjugation(kx, ky, kz);

          T re = val.real();
          T imag = val.imag();

          if(conj) imag=-imag;

          field[index_re] = re;
          field[index_re+1] = imag;

        }

        std::complex<T> getFourierCoefficient(int kx, int ky, int kz)  {
          bool conj;
          size_t index_re;

          std::tie(conj, index_re) = getRealCoeffLocationAndConjugation(kx, ky, kz);

          T re = field[index_re];
          T im = field[index_re+1];
          if(conj)
            im = -im;
          return std::complex<T>(re,im);

        }

        size_t getRequiredDataSize()  {
          // for FFTW3 real<->complex FFTs
          // see http://www.fftw.org/fftw3_doc/Real_002ddata-DFT-Array-Format.html#Real_002ddata-DFT-Array-Format
          return 2*field.getGrid().size2*(field.getGrid().size/2+1);
        }

        void performTransform()  {
          auto &fieldData = field.getDataVector();

          initialise();

          bool transformToFourier = !field.isFourier();

          fftw_plan plan;
          size_t i;

          int res = static_cast<int>(field.getGrid().size);
          double norm = pow(static_cast<double>(res), 1.5);


          if (transformToFourier) {
            padForFFTWRealTransform();
            plan = fftw_plan_dft_r2c_3d(res, res, res,
                                        &fieldData[0],
                                        reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                        FFTW_ESTIMATE);
          } else {
            mirrorFourierModesForFFTWRealTransform();
            plan = fftw_plan_dft_c2r_3d(res, res, res,
                                        reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                        &fieldData[0],
                                        FFTW_ESTIMATE);
          }


          fftw_execute(plan);
          fftw_destroy_plan(plan);

          if(!transformToFourier) {
            unpadAfterFFTWRealTransform();
          }

          using tools::numerics::operator/=;
          fieldData/=norm;



          field.setFourier(!field.isFourier());

        }

      };

      template<>
      class FieldFourierManager<std::complex<double>>  {
      protected:
        using T=double;
        fields::Field<std::complex<T>, T> &field;
        const grids::Grid<T> &grid;
      public:
        FieldFourierManager(fields::Field<std::complex<T>, T> &field) : field(field), grid(field.getGrid()) {

        }

        void setFourierCoefficient(int kx, int ky, int kz, const std::complex<T> &val)  {
          size_t id_k, id_negk;

          id_k = grid.getCellIndex(Coordinate<int>(kx, ky, kz));
          id_negk = grid.getCellIndex(Coordinate<int>(-kx, -ky, -kz));

          field[id_k] = val;
          field[id_negk] = std::conj(val);
        }

        std::complex<T> getFourierCoefficient(int kx, int ky, int kz)  {
          return field[grid.getCellIndex(Coordinate<int>(kx, ky, kz))];
        }

        size_t getRequiredDataSize()  {
          return field.getGrid().size3;
        }

        void performTransform()  {

          auto &fieldData = field.getDataVector();

          initialise();

          fftw_plan plan;
          size_t i;

          int res = static_cast<int>(field.getGrid().size);
          double norm = pow(static_cast<double>(res), 1.5);

          if (!field.isFourier())
            plan = fftw_plan_dft_3d(res, res, res,
                                    reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                    reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                    FFTW_FORWARD, FFTW_ESTIMATE);

          else
            plan = fftw_plan_dft_3d(res, res, res,
                                    reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                    reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                    FFTW_BACKWARD, FFTW_ESTIMATE);


          fftw_execute(plan);
          fftw_destroy_plan(plan);

          using tools::numerics::operator/=;
          fieldData/=norm;

          field.setFourier(!field.isFourier());
        }


      };

      template<typename T>
      int getNyquistModeThatMustBeReal(const grids::Grid<T> &g);

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
      fields::Field<tools::datatypes::ensure_complex<T>, tools::datatypes::strip_complex<T>>
      getComplexFourierField(const fields::Field<T, tools::datatypes::strip_complex<T>> &field) {
        using complex_T = tools::datatypes::ensure_complex<T>;
        using underlying_T = tools::datatypes::strip_complex<T>;
        const grids::Grid<underlying_T> &g = field.getGrid();
        fields::Field<complex_T, underlying_T> out(g);

        std::cerr << "field-size=" << out.getDataVector().size() << std::endl;

        for (size_t i = 0; i < g.size3; ++i) {
          int kx, ky, kz;
          std::tie(kx, ky, kz) = g.getFourierCellCoordinate(i);
          out[i] = field.getFourierCoefficient(kx,ky,kz);
        }
        return out;
      };



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
      int getFourierCellWeight(const fields::Field<std::complex<T>, T> &field, size_t i) {
        return 1;
      }

      template<typename T>
      int getFourierCellWeight(const fields::Field<T, T> &field, size_t i) {
        return 1;
      }


      template<typename T, typename S>
      void performFFT(fields::Field<T, S> &field);

      template<>
      void performFFT(fields::Field<std::complex<double>, double> &field) {

        auto &fieldData = field.getDataVector();

        initialise();

        fftw_plan plan;
        size_t i;

        int res = static_cast<int>(field.getGrid().size);

        double norm = pow(static_cast<double>(res), 1.5);
        size_t len = static_cast<size_t>(res * res);
        len *= res;

        if (!field.isFourier())
          plan = fftw_plan_dft_3d(res, res, res,
                                  reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                  reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                  FFTW_FORWARD, FFTW_ESTIMATE);

        else
          plan = fftw_plan_dft_3d(res, res, res,
                                  reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                  reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                  FFTW_BACKWARD, FFTW_ESTIMATE);


        fftw_execute(plan);
        fftw_destroy_plan(plan);

#pragma omp parallel for schedule(static) private(i)
        for (i = 0; i < len; i++)
          fieldData[i] /= norm;

        field.setFourier(!field.isFourier());


      }


    }
  }
}
#endif // FFTH_INCLUDED
