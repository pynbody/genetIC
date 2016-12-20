//
// Created by Andrew Pontzen on 27/07/15.
//

#ifndef IC_RANDOMFIELDGENERATOR_HPP
#define IC_RANDOMFIELDGENERATOR_HPP

#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "grid.hpp"

template<typename MyFloat>
class RandomFieldGenerator {
protected:
  gsl_rng *randomState;
  const gsl_rng_type *T;
  bool drawInFourierSpace;
  bool reverseRandomDrawOrder;
  MultiLevelField<complex<MyFloat>> & field;


public:
  RandomFieldGenerator(MultiLevelField<complex<MyFloat>> &field_, int seed = 0) :
    field(field_) {
    T = gsl_rng_ranlxs2; // shouldn't this be gsl_rng_ranlxd2 for MyFloat = double? -> it's single precision for compatibility with previous versions!
    randomState = gsl_rng_alloc(T); //this allocates memory for the generator with type T
    gsl_rng_set(randomState, seed);
    drawInFourierSpace = false;
  }

  virtual ~RandomFieldGenerator() {
    gsl_rng_free(randomState);
  }

  using RefFieldType = std::vector<std::complex<MyFloat>> &;
  using FieldType = std::remove_reference_t<RefFieldType>;


  void setDrawInFourierSpace(bool value) {
    drawInFourierSpace = value;
  }

  void setReverseRandomDrawOrder(bool value) {
    reverseRandomDrawOrder = value;
  }

  void seed(int seed) {
    gsl_rng_set(randomState, seed);
  }

  void draw() {
    for(size_t i=0; i<field.getNumLevels(); ++i) {
      auto &fieldOnGrid = field.getFieldOnGrid(i);
      if (drawInFourierSpace) {
        fieldOnGrid.toFourier();
        drawRandomForSpecifiedGridFourier(fieldOnGrid.getGrid(), fieldOnGrid.getDataVector());
      } else {
        drawRandomForSpecifiedGrid(fieldOnGrid.getGrid(), fieldOnGrid.getDataVector());
      }
    }
  }

protected:

  void drawOneFourierMode(const Grid<MyFloat> &g, int k1, int k2, int k3,
                          MyFloat norm, RefFieldType pField_k) {
    size_t id_k, id_negk;

    id_k = g.getCellIndex(Coordinate<int>(k1, k2, k3));
    id_negk = g.getCellIndex(Coordinate<int>(-k1, -k2, -k3));


    // these need to be initialized in explicit order - can't leave it to compilers
    // to choose as they choose differently...
    MyFloat a = norm * gsl_ran_gaussian_ziggurat(randomState, 1.);
    MyFloat b = norm * gsl_ran_gaussian_ziggurat(randomState, 1.);

    if (reverseRandomDrawOrder)
      pField_k[id_k] = std::complex<MyFloat>(b, a);
    else
      pField_k[id_k] = std::complex<MyFloat>(a, b);

    // reality condition:
    pField_k[id_negk] = std::conj(pField_k[id_k]);

  }


  void drawRandomForSpecifiedGrid(const Grid<MyFloat> &g, RefFieldType pField_k) {
    size_t nPartTotal = g.size3;

    MyFloat sigma = sqrt((MyFloat) (nPartTotal));

    cerr << "Drawing random numbers...";

    // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
    for (size_t i = 0; i < nPartTotal; i++) {
      pField_k[i] = gsl_ran_gaussian_ziggurat(randomState, 1.) * sigma;
    }

    // this FFT is kept for historical compatibility, even though it could be made
    // unnecessary by picking random phases directly in fourier space
    fft(pField_k, pField_k, 1);

    pField_k[0] = complex<MyFloat>(0., 0.); // ensure mean==0
    cerr << "done" << endl;
  }

  void drawRandomForSpecifiedGridFourier(const Grid<MyFloat> &g,
                                         RefFieldType pField_k) {

    progress::ProgressBar pb("");
    MyFloat sigma = sqrt((MyFloat) (g.size3));

    cerr << "Drawing random numbers in fourier space..." << endl;
    int ks, k1, k2;

    sigma /= sqrt(2.0);

    // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
    // Do it in square k-shells, in order of increasing |k|, so that
    // resolution can be scaled by factors of 2 and we still get
    // the 'same' field
    for (ks = 0; ks < int(g.size / 2); ks++) {
      pb.setProgress(float(ks*2)/g.size);
      for (k1 = -ks; k1 < ks; k1++) {
        for (k2 = -ks; k2 < ks; k2++) {
          drawOneFourierMode(g, ks, k1, k2, sigma, pField_k);
          drawOneFourierMode(g, k1, ks, k2, sigma, pField_k);
          drawOneFourierMode(g, k1, k2, ks, sigma, pField_k);
        }
      }
    }


  }


};

#endif //IC_RANDOMFIELDGENERATOR_HPP
