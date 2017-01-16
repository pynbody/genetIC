#ifndef _CONSTRAINTAPPLICATOR_HPP
#define _CONSTRAINTAPPLICATOR_HPP

#include "multilevelcontext.hpp"
#include "multilevelfield.hpp"

template<typename DataType, typename T=strip_complex<DataType>>
class ConstraintApplicator {
private:


public:
  std::vector<ConstraintField<DataType>> alphas;
  std::vector<DataType> values;
  std::vector<DataType> existing_values;

  MultiLevelContextInformation<DataType> *underlying;
  OutputField<DataType> *outputField;

  ConstraintApplicator(MultiLevelContextInformation<DataType> *underlying_,
                       OutputField<DataType> *outputField_) : underlying(underlying_), outputField(outputField_) {

  }

  void add_constraint(ConstraintField<DataType> &&alpha, DataType value, DataType existing) {

    alphas.push_back(std::move(alpha));
    values.push_back(value);
    existing_values.push_back(existing);
  }

  void print_covariance() {
    progress::ProgressBar pb("calculating covariance");
    size_t n = alphas.size();
    size_t done = 0;
    std::vector<std::vector<std::complex<T>>> c_matrix(n, std::vector<std::complex<T>>(n, 0));

    // Gram-Schmidt orthogonalization in-place
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j <= i; j++) {
        pb.setProgress( ((float) done * 2) / (n * (1 + n)));
        c_matrix[i][j] = alphas[i].innerProduct(alphas[j]);
        c_matrix[j][i] = c_matrix[i][j];
        done += 1;
      }
    }

    std::cout << std::endl << "cov_matr = [";
    for (size_t i = 0; i < n; i++) {
      std::cout << "[";
      for (size_t j = 0; j < n; j++) {
        std::cout << std::real(c_matrix[i][j]);
        if (j < n - 1) std::cout << ",";
      }
      std::cout << "]";
      if (i < n - 1) std::cout << "," << std::endl;
    }
    std::cout << "]" << std::endl;
  }

  void prepare() {

    size_t n = alphas.size();
    size_t done = 0;

    size_t nCells = underlying->getNumCells();


    std::cout << "v0=[";
    for (size_t i = 0; i < n; i++) {
      std::cout << std::real(existing_values[i]);
      if (i < n - 1)
        std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "v1=[";
    for (size_t i = 0; i < n; i++) {
      std::cout << std::real(values[i]);
      if (i < n - 1)
        std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    // std::cout << "ratio[0] = " << values[0]/existing_values[0] << std::endl;

    // Store transformation matrix, just for display purposes
    std::vector<std::vector<std::complex<T>>> t_matrix(n, std::vector<std::complex<T>>(n, 0));

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        if (i == j) t_matrix[i][j] = 1.0; else t_matrix[i][j] = 0.0;
      }
    }

    // Gram-Schmidt orthogonalization in-place

    progress::ProgressBar pb("orthogonalizing constraints");

    for (size_t i = 0; i < n; i++) {
      auto &alpha_i = alphas[i];
      for (size_t j = 0; j < i; j++) {
        auto &alpha_j = alphas[j];
        pb.setProgress(((float) done * 2) / (n * (1 + n)));
        DataType result = alpha_i.innerProduct(alpha_j);

        alpha_i.addScaled(alpha_j,-result);


        // update display matrix accordingly such that at each step
        // our current values array is equal to t_matrix . input_values
        for (size_t k = 0; k <= j; k++)
          t_matrix[i][k] -= result * t_matrix[j][k];

        // update constraining value
        values[i] -= result * values[j];
        done += 1; // one op for each of the orthognalizing constraints
      }

      // normalize
      DataType norm = sqrt(alpha_i.innerProduct(alpha_i));

      alpha_i/=norm;
      values[i] /= norm;


      for (size_t j = 0; j < n; j++) {
        t_matrix[i][j] /= norm;
      }
      done += 1; // one op for the normalization

    }

    // Now display t_matrix^{dagger} t_matrix, which is the matrix allowing one
    // to form chi^2: Delta chi^2 = d_1^dagger t_matrix^dagger t_matrix d_1 -
    // d_0^dagger t_matrix^dagger t_matrix d_0.



    existing_values.clear();


    for (size_t i = 0; i < n; i++) {
      auto &alpha_i = alphas[i];
      existing_values.push_back(alpha_i.innerProduct(*outputField ));
    }



    //std::cout << "ratio[0] = " << values[0]/existing_values[0] << std::endl;

    std::cout << std::endl << "chi2_matr = [";
    for (size_t i = 0; i < n; i++) {
      std::cout << "[";
      for (size_t j = 0; j < n; j++) {
        std::complex<T> r = 0;
        for (size_t k = 0; k < n; k++) {
          r += std::conj(t_matrix[k][i]) * t_matrix[k][j];
        }
        std::cout << std::real(r) / nCells;
        if (j < n - 1) std::cout << ",";
      }
      std::cout << "]";
      if (i < n - 1) std::cout << "," << std::endl;
    }
    std::cout << "]" << std::endl;


  }

  T get_delta_chi2() {
    T rval = 0.0, v0, v1;
    for (size_t i = 0; i < alphas.size(); i++) {
      v0 = std::abs(existing_values[i]);
      v1 = std::abs(values[i]);
      rval += v1 * v1 - v0 * v0;
    }
    return rval / underlying->getNumCells();
  }


  void applyConstraints() {

    outputField->toFourier();

    for (size_t i = 0; i < alphas.size(); i++) {
      MultiLevelField<DataType> &alpha_i = alphas[i];
      auto dval_i = values[i] - existing_values[i];
      alpha_i.convertToVector();
      alpha_i.toFourier(); // probably already is, but just to be safe
      outputField->addScaled(alpha_i, dval_i);

    }
  }

};


#endif