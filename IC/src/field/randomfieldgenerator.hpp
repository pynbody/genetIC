//
// Created by Andrew Pontzen on 27/07/15.
//

#ifndef IC_RANDOMFIELDGENERATOR_HPP
#define IC_RANDOMFIELDGENERATOR_HPP

#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "src/grid.hpp"

namespace fields {

  template<typename DataType>
  class RandomFieldGenerator {
  protected:
    using FloatType=  strip_complex<DataType>;

    gsl_rng *randomState;
    const gsl_rng_type *randomNumberGeneratorType;
    bool drawInFourierSpace;
    bool reverseRandomDrawOrder;
    bool seeded;
    MultiLevelField <DataType> &field;


  public:
    RandomFieldGenerator(MultiLevelField <DataType> &field_, int seed = 0) :
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

    void seed(int seed) {
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
          drawRandomForSpecifiedGrid(fieldOnGrid.getGrid(), fieldOnGrid.getDataVector());
        }
      }
    }

  protected:


    void drawOneFourierMode(Field<DataType> &field, int k1, int k2, int k3,
                            FloatType norm) {

      // these need to be initialized in explicit order - can't leave it to compilers
      // to choose as they choose differently...
      FloatType a = norm * gsl_ran_gaussian_ziggurat(randomState, 1.);
      FloatType b = norm * gsl_ran_gaussian_ziggurat(randomState, 1.);

      if (reverseRandomDrawOrder)
        fourier::setFourierCoefficient(field, k1, k2, k3, b, a);
      else
        fourier::setFourierCoefficient(field, k1, k2, k3, a, b);

    }


    void drawRandomForSpecifiedGrid(const Grid<FloatType> &g, RefFieldType pField_k) {
      size_t nPartTotal = g.size3;

      FloatType sigma = sqrt((FloatType) (nPartTotal));

      cerr << "Drawing random numbers...";

      // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
      for (size_t i = 0; i < nPartTotal; i++) {
        pField_k[i] = gsl_ran_gaussian_ziggurat(randomState, 1.) * sigma;
      }

      // this FFT is kept for historical compatibility, even though it could be made
      // unnecessary by picking random phases directly in fourier space
      fourier::fft(pField_k, pField_k, 1);

      set_zero(pField_k[0]);

      cerr << "done" << endl;
    }

    void drawRandomForSpecifiedGridFourier(Field<DataType> &field) {
      std::vector<DataType> &vec = field.getDataVector();

      std::fill(vec.begin(), vec.end(), DataType(0));

      const Grid<FloatType> &g = field.getGrid();

      progress::ProgressBar pb("");
      FloatType sigma = sqrt((FloatType) (g.size3));

      cerr << "Drawing random numbers in fourier space..." << endl;
      int ks, k1, k2;

      sigma /= sqrt(2.0);

      // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
      // Do it in square k-shells, in order of increasing |k|, so that
      // resolution can be scaled by factors of 2 and we still get
      // the 'same' field
      for (ks = 0; ks < int(g.size / 2); ks++) {
        pb.setProgress(float(ks * 2) / g.size);
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
