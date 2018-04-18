#ifndef IC_RANDOMFIELDGENERATOR_HPP
#define IC_RANDOMFIELDGENERATOR_HPP

#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "src/simulation/grid/grid.hpp"

namespace fields {

  /** This class handles drawing random white noise across all grids in a multi-level field. */
  template<typename DataType>
  class RandomFieldGenerator {
  protected:
    using FloatType = tools::datatypes::strip_complex<DataType>;

    gsl_rng *randomState;
    const gsl_rng_type *randomNumberGeneratorType;
    bool drawInFourierSpace;
    bool reverseRandomDrawOrder;
    bool seeded;
    MultiLevelField <DataType> &field;


  public:
    RandomFieldGenerator(MultiLevelField <DataType> &field_, unsigned long seed = 0) :
        field(field_) {
      randomNumberGeneratorType = gsl_rng_ranlxs2; // shouldn't this be gsl_rng_ranlxd2 for FloatType = double? -> it's single precision for compatibility with previous versions!
      randomState = gsl_rng_alloc(randomNumberGeneratorType); //this allocates memory for the generator with type T
      gsl_rng_set(randomState, seed);
      drawInFourierSpace = false;
      seeded = false;
    }

    virtual ~RandomFieldGenerator() {
      gsl_rng_free(randomState);
    }

    using RefFieldType = std::vector<DataType> &;
    using FieldType = std::remove_reference_t<RefFieldType>;


    void setDrawInFourierSpace(bool value) {
      drawInFourierSpace = value;
    }

    void setReverseRandomDrawOrder(bool value) {
      reverseRandomDrawOrder = value;
    }

    void seed(unsigned long seed) {
      if (seeded)
        throw std::runtime_error("The random number generator has already been seeded");
      else
        gsl_rng_set(randomState, seed);
      seeded = true;
    }

    void draw() {
      if (!seeded)
        throw std::runtime_error("The random number generator has not been seeded");
      for (size_t i = 0; i < field.getNumLevels(); ++i) {
        auto &fieldOnGrid = field.getFieldForLevel(i);
        if (drawInFourierSpace) {
          fieldOnGrid.toFourier();
          drawRandomForSpecifiedGridFourier(fieldOnGrid);
        } else {
          drawRandomForSpecifiedGrid(fieldOnGrid);
        }
      }
    }

  protected:


    void drawOneFourierMode(Field <DataType> &field, int k1, int k2, int k3,
                            FloatType norm) {

      // these need to be initialized in explicit order - can't leave it to compilers
      // to choose as they choose differently...
      FloatType a = norm * gsl_ran_gaussian_ziggurat(randomState, 1.);
      FloatType b = norm * gsl_ran_gaussian_ziggurat(randomState, 1.);

      if (reverseRandomDrawOrder)
        field.setFourierCoefficient(k1, k2, k3, tools::datatypes::ensure_complex<DataType>(b, a));
      else
        field.setFourierCoefficient(k1, k2, k3, tools::datatypes::ensure_complex<DataType>(a, b));

      int nyquist = int(field.getGrid().size) / 2;

      if (k1 == 0 || k1 == nyquist || k2 == 0 || k2 == nyquist || k3 == 0 || k3 == nyquist) {
        // Due to poor original implementation of drawRandomForSpecifiedGridFourier (which we're now stuck with
        // for historical compatibility), we need to ensure the _last_ mode written to a field (which may, due to
        // symmetries in the Fourier coeffs at 0 or nyquist modes, overwrite a previous draw) persists.
        //
        // The above condition probably catches too many cases, but if there is a risk, let's also explicitly write
        // the related coeff.
        if (reverseRandomDrawOrder)
          field.setFourierCoefficient(-k1, -k2, -k3, tools::datatypes::ensure_complex<DataType>(b, -a));
        else
          field.setFourierCoefficient(-k1, -k2, -k3, tools::datatypes::ensure_complex<DataType>(a, -b));
      }
    }


    void drawRandomForSpecifiedGrid(Field <DataType> &field) {
      /* Draw random white noise in real space.
       *
       *  Kept for historical compatibility, even though the recommended approach is
       * to seed in Fourier space (routine drawRandomForSpecifiedGridFourier, accessed using command
       * seedfourier in the paramfile)
       *
       */

      field.setFourier(false);
      auto &g = field.getGrid();
      auto &fieldData = field.getDataVector();
      size_t nPartTotal = g.size3;

      std::cerr << "Drawing random numbers...";

      // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
      for (size_t i = 0; i < nPartTotal; i++) {
        fieldData[i] = gsl_ran_gaussian_ziggurat(randomState, 1.);
      }

      field.toFourier();
      tools::set_zero(fieldData[0]);

      std::cerr << "done" << std::endl;
    }

    void drawRandomForSpecifiedGridFourier(Field <DataType> &field) {
      /* Draw random white noise in Fourier space. */
      std::vector<DataType> &vec = field.getDataVector();

      std::fill(vec.begin(), vec.end(), DataType(0));

      const grids::Grid<FloatType> &g = field.getGrid();

      tools::progress::ProgressBar pb("");

      std::cerr << "Drawing random numbers in fourier space..." << std::endl;
      int ks, k1, k2;

      FloatType  sigma = 1.0 / sqrt(2.0);

      // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
      // Do it in square k-shells, in order of increasing |k|, so that
      // resolution can be scaled by factors of 2 and we still get
      // the 'same' field
      for (ks = 0; ks < int(g.size / 2); ks++) {
        pb.setProgress(float(ks * ks) * (ks * 8) / g.size3);
        for (k1 = -ks; k1 < ks; k1++) {
          for (k2 = -ks; k2 < ks; k2++) {
            drawOneFourierMode(field, ks, k1, k2, sigma);
            drawOneFourierMode(field, k1, ks, k2, sigma);
            drawOneFourierMode(field, k1, k2, ks, sigma);
          }
        }
      }


    }


  };

}

#endif //IC_RANDOMFIELDGENERATOR_HPP
