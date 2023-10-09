#ifndef FFTH_INCLUDED
#define FFTH_INCLUDED


#include <stdexcept>
#ifdef USE_CUFFT
#include <cufftw.h>
#else
#include <fftw3.h>
#endif

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
    /*! \namespace tools::numerics::fourier
        \brief Provides code that handles fourier transforms, and interfaces with the FFTW library.
    */
    namespace fourier {

      bool fftwThreadsInitialised = false;

      //! Initialises the FFTW threads if they haven't already been initialised
      void initialise() {
        if (fftwThreadsInitialised)
          return;

#ifdef USE_CUFFT
       // nothing to do
       logging::entry() << "Using CUFFT" << std::endl;
#else
#ifdef FFTW_THREADS
        if (fftw_init_threads() == 0)
          throw std::runtime_error("Cannot initialize FFTW threads");
#ifndef _OPENMP
        fftw_plan_with_nthreads(FFTW_THREADS);
  logging::entry() << "Note: " << FFTW_THREADS << " FFTW Threads were initialised" << std::endl;
#else
        int numThreads = omp_get_max_threads();
        bool emitThreadLimitMessage = false;
#if __APPLE__ && __MACH__
#ifndef IGNORE_APPLE_FFTW_THREAD_LIMIT
        if(numThreads>8) {
          emitThreadLimitMessage = true;
          numThreads=8;
        }
#endif
#endif
        fftw_plan_with_nthreads(numThreads);
        if(emitThreadLimitMessage) {
          logging::entry() << std::endl;
          logging::entry()  << "Limiting number of FFTW Threads to " << numThreads << ", because FFTW on Mac OS seems to become slow beyond this point."
            << std::endl;
          logging::entry()
            << "To disable this behaviour, recompile with -DIGNORE_APPLE_FFTW_THREAD_LIMIT" << std::endl;
          logging::entry()
            << "OpenMP parts of the code will still run with " << omp_get_max_threads() << " threads." << std::endl;
          logging::entry() << std::endl;
        } else {
          logging::entry() << "Note: " << numThreads << " FFTW Threads (determined by OpenMP) were initialised"
                    << std::endl;
        }
#endif
#else
        logging::entry() << "Note: FFTW Threads are not enabled" << std::endl;
#endif
#endif
        fftwThreadsInitialised = true;
      }

      /*! \class FieldFourierManagerBase
          \brief Class that handles all operations to do with Fourier transforms used by the code.
      */
      template<typename DataType, typename CoordinateType = tools::datatypes::strip_complex<DataType>>
      class FieldFourierManagerBase {
      protected:
        using ComplexType = tools::datatypes::ensure_complex<DataType>;
        fields::Field<DataType, CoordinateType> &field; //!< Reference to a field to which Fourier operations will be applied.
        const grids::Grid<CoordinateType> &grid; //!< Constant reference to the grid on which the Fourier transform is applied.
        int size; //!< number of points in the Fourier transform
        int nyquistIfEvenElseZero; //!< Zero if odd number of points, otherwise, half the number of points.
        int largestNegativeMode; //!< Largest negative mode for Fourier transforms of this size.


        //! Construct a FieldFourierManagerBase object for a given field
        FieldFourierManagerBase(fields::Field<DataType, CoordinateType> &field) : field(field), grid(field.getGrid()) {
          size = static_cast<int>(grid.size);
          nyquistIfEvenElseZero = size / 2;
          largestNegativeMode = -size / 2;

          if (size % 2 == 0)
            largestNegativeMode += 1;
          else
            nyquistIfEvenElseZero = 0;
        }

        //! Applies the callback function iteratively over a Fourier space field, summing up the result of the function over all cells
        ComplexType
        iterateFourierCellsWithAccumulation(const std::function<ComplexType(int, int, int)> &callback) const {
          field.toFourier();

          CoordinateType global_result_real(0), global_result_imag(0);

#pragma omp parallel for reduction(+:global_result_real, global_result_imag)
          for (int kx = 0; kx < size / 2 + 1; kx++) {
            int ky_lower = largestNegativeMode;
            int cellWeightX = 2;
            if (kx == 0 || kx == nyquistIfEvenElseZero) {
              ky_lower = 0;
              cellWeightX = 1;
            }

            for (int ky = ky_lower; ky < size / 2 + 1; ky++) {
              int cellWeightY = cellWeightX;
              int kz_lower = largestNegativeMode;

              if ((kx == 0 || kx == nyquistIfEvenElseZero) && (ky == 0 || ky == nyquistIfEvenElseZero)) {
                kz_lower = 0;
              }

              if (ky_lower == 0 && ky != 0 && ky != nyquistIfEvenElseZero) {
                cellWeightY *= 2;
              }

              for (int kz = kz_lower; kz < size / 2 + 1; kz++) {
                int cellWeightZ = cellWeightY;
                if (kz_lower == 0 && kz != 0 && kz != nyquistIfEvenElseZero) {
                  cellWeightZ *= 2;
                }
                ComplexType result = callback(kx, ky, kz);
                global_result_real += result.real() * cellWeightZ;
                global_result_imag += result.imag() * cellWeightZ;
              }
            }
          }

          return ComplexType(global_result_real, global_result_imag);

        }

        //! Apply the callback function to all cells (no return type)
        void iterateFourierCells(const std::function<void(int, int, int)> &callback) const {
          // TODO: this is good because we only need to implement the loop over fourier cells once,
          // but bad because functions that are not actually accumulating will waste time accumulating zeros
          auto callback_wrapper = [&callback](int x, int y, int z) -> ComplexType {
            callback(x, y, z);
            return 0;
          };
          iterateFourierCellsWithAccumulation(callback_wrapper);
        }


      public:

        //! Destructor
        virtual ~FieldFourierManagerBase() {

        }

        //! Used to ensure that the Fourier transform of a real quantity will have correctly mirrored Fourier modes
        virtual void ensureFourierModesAreMirrored() {

        }

        /*! \brief Iterate (potentially in parallel) over each Fourier cell, applying the function fn to each cell
            \param fn - The passed function takes arguments (value, kx, ky, kz) where value is the Fourier coeff value
           * at k-mode kx, ky, kz. It returns the new value.

             Updates the corresponding value according to the
           * passed function return value.
           *
           */
        void forEachFourierCell(
          const std::function<ComplexType(ComplexType, CoordinateType, CoordinateType, CoordinateType)> &fn) {

          CoordinateType kMin = grid.getFourierKmin();
          iterateFourierCells([&fn, this, kMin](int kx, int ky, int kz) {
            ComplexType value = getFourierCoefficient(kx, ky, kz);
            value = fn(value, kx * kMin, ky * kMin, kz * kMin);
            setFourierCoefficient(kx, ky, kz, value);
          });
        }


        /*! \brief Iterate (potentially in parallel) over each Fourier cell, applying function fn which returns no result
            \param fn - The passed function takes arguments (value, kx, ky, kz) where value is the Fourier coeff value
           * at k-mode kx, ky, kz.
           */
        void forEachFourierCell(
          const std::function<void(ComplexType, CoordinateType, CoordinateType, CoordinateType)> &fn) const {
          CoordinateType kMin = grid.getFourierKmin();
          iterateFourierCells([&fn, this, kMin](int kx, int ky, int kz) {
            ComplexType value = getFourierCoefficient(kx, ky, kz);
            fn(value, kx * kMin, ky * kMin, kz * kMin);
          });
        }

        /*! \brief Iterate (potentially in parallel) over each Fourier cell.
            \param fn - The passed function takes arguments (value, kx, ky, kz) where value is the Fourier coeff value
           * at k-mode corresponding to the integer wavenumbers kx*grid.getFourierKmin() etc
           */
        void forEachFourierCellInt(const std::function<void(ComplexType, int, int, int)> &fn) const {
          CoordinateType kMin = grid.getFourierKmin();
          iterateFourierCells([&fn, this, kMin](int kx, int ky, int kz) {
            ComplexType value = getFourierCoefficient(kx, ky, kz);
            fn(value, kx, ky, kz);
          });
        }

        /*! \brief Iterate over, and update the value in, each Fourier cell. Possibly executed in parallel.
            \param fn - The passed function takes arguments (value, kx, ky, kz) where value is the Fourier coeff value
           * at k-mode corresponding to the integer wavenumbers kx*grid.getFourierKmin() etc.
           * It must return the updated value.
           */
        void forEachFourierCellInt(const std::function<ComplexType(ComplexType, int, int, int)> &fn) {
          CoordinateType kMin = grid.getFourierKmin();
          iterateFourierCells([&fn, this, kMin](int kx, int ky, int kz) {
            ComplexType value = getFourierCoefficient(kx, ky, kz);
            auto retval = fn(value, kx, ky, kz);
            setFourierCoefficient(kx, ky, kz, retval);
          });
        }

        /*! \brief Iterate (potentially in parallel) and accumulate a complex number over each Fourier cell.
            \param callback - The passed function takes arguments (value, kx, ky, kz). The return value is accumulated.
           */
        ComplexType
        accumulateForEachFourierCell(const std::function<ComplexType(ComplexType, int, int, int)> &callback) const {
          field.toFourier();
          auto result = iterateFourierCellsWithAccumulation([&callback, this](int kx, int ky, int kz) -> ComplexType {
            ComplexType value = getFourierCoefficient(kx, ky, kz);
            return callback(value, kx, ky, kz);
          });

          return result;

        }


        //! Takes a function outputting a tuple of three complex numbers, and iterates it over the Fourier grid, to create three Fourier fields
        auto generateNewFourierFields(
          const std::function<std::tuple<ComplexType, ComplexType, ComplexType>(ComplexType, CoordinateType,
                                                                                CoordinateType,
                                                                                CoordinateType)> &fn) {
          using Field = fields::Field<DataType, CoordinateType>;
          // TODO: ought to be possible to generalise away from ugly 3D-specific case using template programming
          auto ret1 = std::make_shared<Field>(grid);
          auto ret2 = std::make_shared<Field>(grid);
          auto ret3 = std::make_shared<Field>(grid);
          CoordinateType kMin = grid.getFourierKmin();
          iterateFourierCells([&fn, this, kMin, ret1, ret2, ret3](int kx, int ky, int kz) {
            ComplexType v1, v2, v3;
            ComplexType value = getFourierCoefficient(kx, ky, kz);
            std::tie(v1, v2, v3) = fn(value, kx * kMin, ky * kMin, kz * kMin);
            ret1->setFourierCoefficient(kx, ky, kz, v1);
            ret2->setFourierCoefficient(kx, ky, kz, v2);
            ret3->setFourierCoefficient(kx, ky, kz, v3);
          });

          return std::make_tuple(ret1, ret2, ret3);

        }


      public:

        //! Sets the specified Fourier coefficient to the specified value
        virtual void setFourierCoefficient(int kx, int ky, int kz, const ComplexType &val) = 0;

        //! Returns the value of the Fourier coefficient at the specified integer wave-numbers
        virtual ComplexType getFourierCoefficient(int kx, int ky, int kz) const = 0;
      };

      /*! \class FieldFourierManager<double>
          \brief Fourier manager class, specialised for real fields
      */
      template<typename T>
      class FieldFourierManager<T, T> : public FieldFourierManagerBase<T, T> {
      protected:
        int size; //!< Number of elements in the set to apply discrete Fourier transform to.
        size_t compressed_size; //!< Compressed size, exploiting symmetry of real discrete Fourier transforms.
        fftw_plan forwardPlan; //!< Method used for going from real space to Fourier space, if T is double
        fftw_plan reversePlan; //!< Method used for going from Fourier space to real space, if T is double
        fftwf_plan forwardPlanFloat; //!< Method used for going from real space to Fourier space, if T is float
        fftwf_plan reversePlanFloat; //!< Method used for going from Fourier space to real space, if T is float

        //! Re-organises the wave-numbers to lie in the positive quadrant, and returns to a linear index (and whether we conjugated the field)
        auto getRealCoeffLocationAndConjugation(int kx, int ky, int kz) const {
          bool conjugate = false;
          if (kz < 0) {
            conjugate = true;
            kx = -kx;
            ky = -ky;
            kz = -kz;
          }
          if (ky < 0) {
            ky = size + ky;
          }
          if (kx < 0) {
            kx = size + kx;
          }
          size_t logical_index = kz + compressed_size * ky + compressed_size * size_t(size * kx);
          size_t index_re = 2 * logical_index;
          assert(index_re + 1 < this->field.getDataVector().size());

          return std::make_tuple(conjugate, index_re);
        }

      public:
        /*! \brief Ensures mirrored Fourier modes for real fields.

            For conceptual ease, the FFTW real transformations contain some duplicate data
            see  http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html

            we always address positive fourier modes where there is an ambiguity, but it is
            not clear that FFTW does the same, so before a transform back to real space we need
            to fill in the missing modes.
        */
        void ensureFourierModesAreMirrored() override {


          bool unused;

#pragma omp parallel for
          for (int kx = 0; kx < size / 2 + 1; ++kx) {
            size_t loc_source;
            size_t loc_dest;

            // N.B. take y loop in reverse order so that it copes with zero and nyquist kx modes naturally
            // (where we MUST copy +ve modes into -ve modes, not the opposite)
            for (int ky = size / 2; ky > -size / 2; --ky) {
              std::tie(unused, loc_source) = getRealCoeffLocationAndConjugation(kx, ky, 0);
              std::tie(unused, loc_dest) = getRealCoeffLocationAndConjugation(-kx, -ky, 0);
              this->field[loc_dest] = this->field[loc_source];
              this->field[loc_dest + 1] = -this->field[loc_source + 1];

              // on an odd-sized grid, the following is a null op. On an even sized-grid, it sorts out the kz
              // nyquist mode.
              std::tie(unused, loc_source) = getRealCoeffLocationAndConjugation(kx, ky, this->nyquistIfEvenElseZero);
              std::tie(unused, loc_dest) = getRealCoeffLocationAndConjugation(-kx, -ky, this->nyquistIfEvenElseZero);
              this->field[loc_dest] = this->field[loc_source];
              this->field[loc_dest + 1] = -this->field[loc_source + 1];
            }
          }
        }

      protected:

        /*! \brief Convert to the FFTW real format

          More precisely, converts from our internal array format (which is just a standard column-major rep)
          to the FFTW real-packed version which has padding every grid_size

          see http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html
        */
        void padForFFTWRealTransform() {

          size_t source_range_end = this->grid.size3;
          size_t padding_amount = compressed_size * 2 - this->grid.size;
          size_t target_range_end =
            this->grid.size * this->grid.size * 2 * compressed_size -
            padding_amount;
          auto &data = this->field.getDataVector();

          while (source_range_end > this->grid.size) {
            size_t target_range_start = target_range_end - this->grid.size;
            size_t source_range_start = source_range_end - this->grid.size;
            std::copy_backward(&data[source_range_start], &data[source_range_end], &data[target_range_end]);
            target_range_end = target_range_start - padding_amount;
            source_range_end = source_range_start;
          }


        }

        //! Reverse transformation of padForFFTWRealTransform
        void unpadAfterFFTWRealTransform() {

          size_t source_range_start = 0;
          size_t target_range_start = 0;
          size_t source_max = getRequiredDataSize();

          size_t padding_amount = compressed_size * 2 - this->grid.size;
          auto &data = this->field.getDataVector();

          while (source_range_start < source_max) {
            size_t source_range_end = source_range_start + this->grid.size;
            std::copy(&data[source_range_start], &data[source_range_end], &data[target_range_start]);
            source_range_start = source_range_end + padding_amount;
            target_range_start += this->grid.size;
          }

        }

      public:
        //! Constructor from a real field
        FieldFourierManager(fields::Field<T, T> &field) : FieldFourierManagerBase<T, T>(field) {
          size = static_cast<int>(this->grid.size);
          compressed_size = this->grid.size / 2 + 1;
          forwardPlan = nullptr;
          reversePlan = nullptr;
          forwardPlanFloat = nullptr;
          reversePlanFloat = nullptr;
        }

        //! Destructor
        virtual ~FieldFourierManager() {
          if (forwardPlan != nullptr) {
            fftw_destroy_plan(forwardPlan);
            forwardPlan = nullptr;
          }
          if (reversePlan != nullptr) {
            fftw_destroy_plan(reversePlan);
            reversePlan = nullptr;
          }
          if (forwardPlanFloat != nullptr) {
            fftwf_destroy_plan(forwardPlanFloat);
            forwardPlanFloat = nullptr;
          }
          if (reversePlanFloat != nullptr) {
            fftwf_destroy_plan(reversePlanFloat);
            reversePlanFloat = nullptr;
          }
        }

        //! Sets the specified Fourier coefficient to val (accounting for mirrored Fourier modes as real field)
        void setFourierCoefficient(int kx, int ky, int kz, const std::complex<T> &val) override {
          bool conj;
          size_t index_re;

          std::tie(conj, index_re) = getRealCoeffLocationAndConjugation(kx, ky, kz);

          T re = val.real();
          T imag = val.imag();

          if (conj) imag = -imag;

          this->field[index_re] = re;
          this->field[index_re + 1] = imag;

        }

        //! Returns the specified Fourier coefficient (accounting for mirrored Fourier modes as real field)
        std::complex<T> getFourierCoefficient(int kx, int ky, int kz) const override {
          bool conj;
          size_t index_re;

          std::tie(conj, index_re) = getRealCoeffLocationAndConjugation(kx, ky, kz);

          T re = this->field[index_re];
          T im = this->field[index_re + 1];
          if (conj)
            im = -im;
          return std::complex<T>(re, im);

        }

        //! Returns the size of the data storage needed to store a real Fourier transform (less than for a generic FT)
        size_t getRequiredDataSize() {
          // for FFTW3 real<->complex FFTs
          // see http://www.fftw.org/fftw3_doc/Real_002ddata-DFT-Array-Format.html#Real_002ddata-DFT-Array-Format
          return 2 * this->field.getGrid().size2 * (
            this->field.getGrid().size / 2 + 1);
        }

        //! Performs the Fourier transform, interfacing with FFTW
        void performTransform() {
          auto &fieldData = this->field.getDataVector();

          initialise();

          bool transformToFourier = !this->field.isFourier();

          fftw_plan plan = nullptr;
          fftwf_plan planFloat = nullptr; // only one of these two will be used

          int res = static_cast<int>(this->field.getGrid().size);
          T norm = pow(static_cast<T>(res), 1.5);

          std::tie(plan, planFloat) = getFFTWPlan(transformToFourier);

          if(transformToFourier) {
            padForFFTWRealTransform();
          } else {
            ensureFourierModesAreMirrored();
          }

          if(plan!=nullptr)
            fftw_execute(plan);
          else if(planFloat!=nullptr)
            fftwf_execute(planFloat);
          else
            throw std::runtime_error("No plan available for FFTW transform");


          if (!transformToFourier) {
            unpadAfterFFTWRealTransform();
          }

          using tools::numerics::operator/=;
          fieldData /= norm;


          this->field.setFourier(!this->field.isFourier());

        }

      private:
        auto getFFTWPlan(bool transformToFourier) {
          auto &fieldData = this->field.getDataVector();
          int res = static_cast<int>(this->field.getGrid().size);
          fftw_plan plan(nullptr);
          fftwf_plan planFloat(nullptr);

          if (transformToFourier) {
            if (forwardPlan == nullptr && forwardPlanFloat == nullptr) {
              if(std::is_same<T,double>::value) {
                forwardPlan = fftw_plan_dft_r2c_3d(res, res, res, reinterpret_cast<double*>(&fieldData[0]),
                                                   reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                                   FFTW_ESTIMATE);
              } else {
                forwardPlanFloat = fftwf_plan_dft_r2c_3d(res, res, res, reinterpret_cast<float*>(&fieldData[0]),
                                                         reinterpret_cast<fftwf_complex *>(&fieldData[0]),
                                                         FFTW_ESTIMATE);
              }
            }

            // one, but only one, of these will be nullptr
            plan = forwardPlan;
            planFloat = forwardPlanFloat;

          } else {
            if(reversePlan == nullptr && reversePlanFloat == nullptr) {
              if(std::is_same<T,double>::value) {
                reversePlan = fftw_plan_dft_c2r_3d(res, res, res,
                                                   reinterpret_cast<fftw_complex *>(&fieldData[0]),
                                                   reinterpret_cast<double*>(&fieldData[0]),
                                                   FFTW_ESTIMATE);
              } else {
                reversePlanFloat = fftwf_plan_dft_c2r_3d(res, res, res,
                                                         reinterpret_cast<fftwf_complex *>(&fieldData[0]),
                                                         reinterpret_cast<float*>(&fieldData[0]),
                                                         FFTW_ESTIMATE);
              }
            }
            // one, but only one, of these will be nullptr
            plan = reversePlan;
            planFloat = reversePlanFloat;
          }
          return std::make_pair(plan, planFloat);
        }


      };

      //! FieldFourierManager specialisation to deal with Fourier transforms of complex fields.
      template<typename T>
      class FieldFourierManager<std::complex<T>, T> : public FieldFourierManagerBase<std::complex<T>, T> {
      public:
        //! Constructor from a complex field
        FieldFourierManager(fields::Field<std::complex<T>, T> &field) : FieldFourierManagerBase<std::complex<T>, T>(field) {

        }

        //! Sets the specified Fourier coefficient to val
        void setFourierCoefficient(int kx, int ky, int kz, const std::complex<T> &val) {
          size_t id_k, id_negk;

          id_k = this->grid.getIndexFromCoordinate(Coordinate<int>(kx, ky, kz));
          id_negk = this->grid.getIndexFromCoordinate(Coordinate<int>(-kx, -ky, -kz));

          this->field[id_k] = val;
          this->field[id_negk] = std::conj(val);
        }

        //! Returns the specified Fourier coefficient
        std::complex<T> getFourierCoefficient(int kx, int ky, int kz) const {
          return this->field[this->grid.getIndexFromCoordinate(Coordinate<int>(kx, ky, kz))];
        }

        //! Returns space required to store Fourier information (always the size of the full grid for Fourier transforms of complex fields)
        size_t getRequiredDataSize() {
          return this->field.getGrid().size3;
        }

        //! Performs the Fourier transform operation, in this cases assuming the field is generic (ie, complex)
        void performTransform() {

          auto &fieldData = this->field.getDataVector();

          initialise();

          fftw_plan plan;

          int res = static_cast<int>(this->field.getGrid().size);
          double norm = pow(static_cast<double>(res), 1.5);

          if (!this->field.isFourier())
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
          fieldData /= norm;

          this->field.setFourier(!this->field.isFourier());
        }


      };

      //! A dummy specialisation expressing that fourier transforms over boolean masks can't be done!
      template<typename T>
      class FieldFourierManager<char, T> : public FieldFourierManagerBase<char, T> {
        using ComplexType = std::complex<char>;
      public:

        void ensureFourierModesAreMirrored() override {

        }

        FieldFourierManager(fields::Field<char, T> &field) : FieldFourierManagerBase<char, T>(field) {

        }

        size_t getRequiredDataSize() {
          return this->field.getGrid().size3;
        }

        void performTransform() {
          throw std::runtime_error("Boolean fields cannot be Fourier transformed");

        }

        void setFourierCoefficient(int kx, int ky, int kz, const ComplexType &val) override {
          throw std::runtime_error("Boolean fields cannot be Fourier transformed");
        }

        ComplexType getFourierCoefficient(int kx, int ky, int kz) const override {
          throw std::runtime_error("Boolean fields cannot be Fourier transformed");
        }

      };

      //! Returns half the number of elements in the grid if even, and an arbitrary large number otherwise.
      template<typename T>
      int getNyquistModeThatMustBeReal(const grids::Grid<T> &g) {
        int even_nyquist;
        if (g.size % 2 == 0)
          even_nyquist = int(g.size) / 2;
        else
          even_nyquist = int(g.size) * 100; // arbitrary large number that will not be seen

        return even_nyquist;
      }


      //! Creates a copy of the field associated to the Fourier manager
      template<typename T>
      fields::Field<tools::datatypes::ensure_complex<T>, tools::datatypes::strip_complex<T>>
      getComplexFourierField(const fields::Field<T, tools::datatypes::strip_complex<T>> &field) {
        using complex_T = tools::datatypes::ensure_complex<T>;
        using underlying_T = tools::datatypes::strip_complex<T>;
        const grids::Grid<underlying_T> &g = field.getGrid();
        fields::Field<complex_T, underlying_T> out(g);

        field.forEachFourierCellInt([&out](tools::datatypes::ensure_complex<T> value, int kx, int ky, int kz) {
          out.setFourierCoefficient(kx, ky, kz, value);
        });

        return out; // TODO - this seems quite inefficient, because it returns the field by value (after already creating a copy of it once already!)
      };

    }
  }
}
#endif // FFTH_INCLUDED
