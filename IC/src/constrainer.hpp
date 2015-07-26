
#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <cassert>
#include <vector>

//#include <Eigen/Dense>


template<typename T>
class UnderlyingField {
private:
    std::vector<std::shared_ptr<Grid<T>>> pGrid;
    std::vector<std::complex<T> *> C0s;
    std::vector<T> weights;

protected:
    std::vector<size_t> Ns;
    std::vector<size_t> cumu_Ns;
    size_t Ntot;
    size_t n_components;

    void id_to_component_id(size_t i, size_t &component, size_t &component_i) {
        component=0;
        while(component<n_components && i>Ns[component]) {
            i-=Ns[component];
            component++;
        }
        if(i>Ns[component])
            throw std::runtime_error("ID out of range when mapping into underlying fields in id_to_component_id");
        component_i = i;
    }

    UnderlyingField(size_t N) {
        n_components=1;
        Ns.push_back(N);
        Ntot = N;
    }

public:
    UnderlyingField(std::complex<T> *C0,  std::shared_ptr<Grid<T>> pG, size_t N) {
        n_components=0;
        Ntot=0;
        add_component(C0,pG,N);
    }

    virtual ~UnderlyingField() { }

    void add_component(std::complex<T> *C0,  std::shared_ptr<Grid<T>> pG, size_t N) {
        C0s.push_back(C0);
        if(pGrid.size()==0) {
            weights.push_back(1.0);
        } else {
            weights.push_back(pow(pG->dx/pGrid[0]->dx,3.0));
        }
        pGrid.push_back(pG);
        Ns.push_back(N);
        cumu_Ns.push_back(Ntot);
        Ntot+=N;
        n_components+=1;
    }

    long get_ntot() const {
        return Ntot;
    }


    //
    // SPECIFIC calculations for the underlying field
    //

    virtual std::complex<T> cov(const std::vector<std::complex<T>> &vec, size_t i) {
        // returns Sum_j C[i,j] vec[j]
        // Here, though, C0 is diagonal
        size_t comp;
        size_t j;
        id_to_component_id(i,comp,j);
        return vec[i]*C0s[comp][j]*weights[comp];
    }

    virtual std::complex<T> cov(const std::vector<std::complex<T>> &vec, size_t comp, size_t j) {
        // returns Sum_i C[j,i] vec[i]
        // Here, though, C0 is diagonal
        size_t overall_j = j+cumu_Ns[comp];
        return C0s[comp][j]*vec[overall_j]*weights[comp];
    }

    //
    // GENERAL calculations that can run on constrained fields too
    //

    virtual std::complex<T> v1_dot_y(const std::vector<std::complex<T>> &v1, bool fourier=true){
        // Calculate v1^{\dagger} y where y is the underlying realization
        T res_real(0), res_imag(0);


        for(size_t comp=0; comp<n_components; comp++) {
            const auto & field = fourier?pGrid[comp]->getFieldFourier():pGrid[comp]->getFieldReal();

            T weight = weights[comp];

            #pragma omp parallel for reduction(+:res_real,res_imag)
            for(size_t i=0; i<Ns[comp]; i++) {
                std::complex<T> res = conj(v1[i+cumu_Ns[comp]])*field[i]*weight;
                // accumulate separately - OMP doesn't support complex number reduction :-(
                res_real+=std::real(res);
                res_imag+=std::imag(res);
            }
        }
        return std::complex<T>(res_real,res_imag);
    }

    virtual std::complex<T> v1_dot_v2(const std::vector<std::complex<T>> &v1,
                                      const std::vector<std::complex<T>> &v2){
        // Calculate v1^{\dagger} v2
        T res_real(0), res_imag(0);

        for(size_t comp=0; comp<n_components; comp++) {
            T weight = weights[comp];
            // #pragma omp parallel for reduction(+:res_real,res_imag)
            for(size_t i=0; i<Ns[comp]; i++) {

                std::complex<T> res = conj(v1[i+cumu_Ns[comp]])*v2[i+cumu_Ns[comp]]*weight;
                // accumulate separately - OMP doesn't support complex number reduction :-(
                res_real+=std::real(res);
                res_imag+=std::imag(res);
            }
        }
        return std::complex<T>(res_real,res_imag);
    }

    virtual std::complex<T> v1_cov_v2(const std::vector<std::complex<T>> &v1,
                                      const std::vector<std::complex<T>> &v2){
        // Calculate v1^{\dagger} C v2
        T res_real(0), res_imag(0);
        for(size_t comp=0; comp<n_components; comp++) {
            T weight = weights[comp];
            #pragma omp parallel for reduction(+:res_real,res_imag)
            for(size_t i=0; i<Ns[comp]; i++) {
                std::complex<T> res = conj(v1[i+cumu_Ns[comp]])*cov(v2,comp,i)*weight;
                // accumulate separately - OMP doesn't support complex number reduction :-(
                res_real+=std::real(res);
                res_imag+=std::imag(res);
            }
        }
        return std::complex<T>(res_real,res_imag);
    }

    virtual std::complex<T> v_cov_v(const std::vector<std::complex<T>> &v){
        return v1_cov_v2(v,v);
    }


    virtual void get_realization(std::complex<T> *r) {
        // Get the realization and store it in the specified array
        size_t j=0;
        for(size_t comp=0; comp<n_components; comp++) {
            auto field = pGrid[comp]->getFieldFourier();

            for(size_t i=0; i<Ns[comp]; i++) {
                r[j]=field[i];
                j++;
            }
        }
        assert(j==Ntot);
    }

    virtual void get_mean(std::complex<T> *r) {
        // Get the mean of all realizations and store it in the specified array
        for(size_t i=0; i<Ntot; i++) {
            if(i%1000==0) progress("get mean x", float(i)/Ntot);
            r[i]=0;
        }
        end_progress();
    }

    virtual T get_delta_chi2() {
        return 0;
    }

};

template<typename T>
class MultiConstrainedField: public UnderlyingField<T>
{
private:


public:
    std::vector<std::vector<std::complex<T>>> alphas;
    std::vector<std::complex<T> > values;
    std::vector<std::complex<T> > existing_values;

    UnderlyingField<T>* underlying;

    MultiConstrainedField(UnderlyingField<T>* underlying_) : UnderlyingField<T>(underlying_->get_ntot()),
    underlying(underlying_)
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
                std::complex<T> result = underlying->v1_cov_v2(alpha_j,alpha_i);

                #pragma omp parallel for
                for(size_t k=0; k<this->Ntot; k++) {
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
            std::complex<T> norm = sqrt(underlying->v_cov_v(alpha_i));

            #pragma omp parallel for
            for(size_t k=0; k<this->Ntot; k++) {
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
                std::cout << std::real(r)/this->Ntot;
                if(j<n-1) std::cout << ",";
            }
            std::cout << "]";
            if(i<n-1) std::cout << "," << std::endl;
        }
        std::cout << "]" << std::endl;



    }

    T get_delta_chi2() override {
        T rval=0.0, v0, v1;
        for(size_t i=0; i<alphas.size(); i++) {
            v0 = std::abs(existing_values[i]);
            v1 = std::abs(values[i]);
            rval+=v1*v1-v0*v0;
        }
        return rval/this->Ntot;
    }

    std::complex<T> v1_dot_y(const std::vector<std::complex<T>> &v1, bool fourier=true) override {
        return underlying->v1_dot_y(v1, fourier);
    }

    std::complex<T> v1_cov_v2(const std::vector<std::complex<T>> &v1, const std::vector<std::complex<T>> &v2) override {
        return underlying->v1_cov_v2(v1,v2);
    }

    void get_realization(std::complex<T> *r) override {
        size_t n=alphas.size();
        underlying->get_realization(r);



        for(size_t i=0; i<n; i++) {
            auto & alpha_i=alphas[i];
            std::complex<T> dval_i = values[i]-existing_values[i];
            progress("get constrained y", float(i)/n);

            #pragma omp parallel for
            for(size_t j=0; j<this->Ntot; j++) {
                r[j]+=dval_i*underlying->cov(alpha_i,j);
            }
        }

        end_progress();
    }
};
