#ifndef _CONSTRAINTAPPLICATOR_HPP
#define _CONSTRAINTAPPLICATOR_HPP

#include "multilevelfieldmanager.hpp"

template<typename T>
class ConstraintApplicator
{
private:


public:
    std::vector<std::vector<std::complex<T>>> alphas;
    std::vector<std::complex<T> > values;
    std::vector<std::complex<T> > existing_values;

    MultiLevelFieldManager<T>* underlying;

    ConstraintApplicator(MultiLevelFieldManager<T>* underlying_) : underlying(underlying_)
    {

    }

    void add_constraint(std::vector<std::complex<T>> && alpha, std::complex<T> value, std::complex<T> existing) {

        alphas.push_back(std::move(alpha));
        values.push_back(value);
        existing_values.push_back(existing);
    }

    void print_covariance() {
        size_t n=alphas.size();
        size_t done=0;
        std::vector<std::vector<std::complex<T>>> c_matrix(n,std::vector<std::complex<T>>(n,0));

        // Gram-Schmidt orthogonalization in-place
        for(size_t i=0; i<n; i++) {
            auto & alpha_i=alphas[i];
            for(size_t j=0; j<=i; j++) {
                auto & alpha_j=alphas[j];
                progress("calculating covariance", ((float)done*2)/(n*(1+n)));
                c_matrix[i][j]=underlying->v1_cov_v2(alpha_j,alpha_i);
                c_matrix[j][i]=c_matrix[i][j];
                done+=1;
            }
        }
        end_progress();

        std::cout << std::endl << "cov_matr = [";
        for(size_t i=0; i<n; i++) {
            std::cout << "[";
            for(size_t j=0; j<n; j++) {
                std::cout << std::real(c_matrix[i][j]);
                if(j<n-1) std::cout << ",";
            }
            std::cout << "]";
            if(i<n-1) std::cout << "," << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    void prepare() {

        size_t n=alphas.size();
        size_t done=0;

        size_t nCells = underlying->getNumCells();


        std::cout << "v0=[";
        for(size_t i=0; i<n; i++) {
            std::cout << std::real(existing_values[i]);
            if(i<n-1)
                std::cout << ", ";
        }
        std::cout << "]"<<std::endl;

        std::cout << "v1=[";
        for(size_t i=0; i<n; i++) {
            std::cout << std::real(values[i]);
            if(i<n-1)
                std::cout << ", ";
        }
        std::cout << "]"<<std::endl;

        // std::cout << "ratio[0] = " << values[0]/existing_values[0] << std::endl;

        // Store transformation matrix, just for display purposes
        std::vector<std::vector<std::complex<T>>> t_matrix(n,std::vector<std::complex<T>>(n,0));

        for(size_t i=0; i<n; i++) {
            for(size_t j=0; j<n; j++) {
                if(i==j) t_matrix[i][j]=1.0; else t_matrix[i][j]=0.0;
            }
        }

        // Gram-Schmidt orthogonalization in-place
        for(size_t i=0; i<n; i++) {
            auto & alpha_i=alphas[i];
            for(size_t j=0; j<i; j++) {
                auto & alpha_j=alphas[j];
                progress("orthogonalizing constraints", ((float)done*2)/(n*(1+n)));
                std::complex<T> result = underlying->v1_cov_v2_with_pseudo_crossterm(alpha_j,alpha_i);

#pragma omp parallel for
                for(size_t k=0; k<nCells; k++) {
                    alpha_i[k]-=result*alpha_j[k];
                }

                // update display matrix accordingly such that at each step
                // our current values array is equal to t_matrix . input_values
                for(size_t k=0; k<=j; k++)
                    t_matrix[i][k]-=result*t_matrix[j][k];

                // update constraining value
                values[i]-=result*values[j];
                done+=1; // one op for each of the orthognalizing constraints
            }

            // normalize
            std::complex<T> norm = sqrt(underlying->v_cov_v_with_pseudo_crossterm(alpha_i));

#pragma omp parallel for
            for(size_t k=0; k<nCells; k++) {
                alpha_i[k]/=norm;
            }
            values[i]/=norm;



            for(size_t j=0;j<n; j++) {
                t_matrix[i][j]/=norm;
            }
            done+=1; // one op for the normalization

        }
        end_progress();

        // Now display t_matrix^{dagger} t_matrix, which is the matrix allowing one
        // to form chi^2: Delta chi^2 = d_1^dagger t_matrix^dagger t_matrix d_1 -
        // d_0^dagger t_matrix^dagger t_matrix d_0.



        existing_values.clear();
        for(size_t i=0; i<n; i++) {
            auto & alpha_i=alphas[i];
            progress("calculating existing means",((float)i)/n);
            existing_values.push_back(underlying->v1_dot_y(alpha_i));
        }
        end_progress();



        //std::cout << "ratio[0] = " << values[0]/existing_values[0] << std::endl;

        std::cout << std::endl << "chi2_matr = [";
        for(size_t i=0; i<n; i++) {
            std::cout << "[";
            for(size_t j=0; j<n; j++) {
                std::complex<T> r=0;
                for(size_t k=0;k<n; k++) {
                    r+=std::conj(t_matrix[k][i])*t_matrix[k][j];
                }
                std::cout << std::real(r)/nCells;
                if(j<n-1) std::cout << ",";
            }
            std::cout << "]";
            if(i<n-1) std::cout << "," << std::endl;
        }
        std::cout << "]" << std::endl;



    }

    T get_delta_chi2()  {
        T rval=0.0, v0, v1;
        for(size_t i=0; i<alphas.size(); i++) {
            v0 = std::abs(existing_values[i]);
            v1 = std::abs(values[i]);
            rval+=v1*v1-v0*v0;
        }
        return rval/underlying->getNumCells();
    }


    void applyConstraints() {

        for(size_t i=0; i<alphas.size(); i++) {
            auto & alpha_i = alphas[i];
            auto dval_i = values[i] - existing_values[i];

            underlying->forEachCellOfEachLevel([dval_i, &alpha_i,this](size_t level, size_t j,
                                                                  size_t overall_field_j,
                                                                  std::vector<std::complex<T>> &field) {
                field[j]+=dval_i*underlying->cov(alpha_i,overall_field_j);

            });
        }
    }

    void get_realization(std::complex<T> *r) {
        size_t nCells = underlying->getNumCells();
        size_t n=alphas.size();
        underlying->get_realization(r);

        for(size_t i=0; i<n; i++) {
            auto & alpha_i=alphas[i];
            std::complex<T> dval_i = values[i]-existing_values[i];
            progress("get constrained y", float(i)/n);


#pragma omp parallel for
            for(size_t j=0; j<nCells; j++) {
                r[j]+=dval_i*underlying->cov(alpha_i,j);
            }
        }

        end_progress();

    }
};


#endif