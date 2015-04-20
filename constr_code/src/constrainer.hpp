
#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <cassert>
//#include <Eigen/Dense>

void progress(const char* message, float progress) {
    int barWidth = 70;

    std::cout << message << " [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %         \r";
    std::cout.flush();
}

void end_progress() {
    std::cout << "                                                                                                                " << std::endl;
}

template<typename T>
class UnderlyingField {
private:
    std::vector<std::complex<T> *> y0s;
    std::vector<std::complex<T> *> C0s;

protected:
    std::vector<long int> Ns;
    std::vector<long int> cumu_Ns;
    long int Ntot;
    int n_components;

    void id_to_component_id(long int i, int &component, long int &component_i) {
        component=0;
        while(component<n_components && i-Ns[component]>0) {
            i-=Ns[component];
            component++;
        }
        if(i>Ns[component])
            throw std::runtime_error("ID out of range when mapping into underlying fields in id_to_component_id");
        component_i = i;
    }

    UnderlyingField(long int N) {
        n_components=1;
        Ns.push_back(N);
        Ntot = N;
    }

public:
    UnderlyingField(std::complex<T> *C0,  std::complex<T> *y0, long int N) {
        n_components=0;
        Ntot=0;
        add_component(C0,y0,N);
    }

    void add_component(std::complex<T> *C0,  std::complex<T> *y0, long int N) {
        C0s.push_back(C0);
        y0s.push_back(y0);
        Ns.push_back(N);
        cumu_Ns.push_back(Ntot);
        Ntot+=N;
        n_components+=1;
    }


    //
    // SPECIFIC calculations for the underlying field
    //

    virtual std::complex<T> cov(std::complex<T> *vec, long int i) {
        // returns Sum_j C[i,j] vec[j]
        // Here, though, C0 is diagonal
        int comp;
        long j;
        id_to_component_id(i,comp,j);
        return C0s[comp][j]*vec[i];
    }

    virtual std::complex<T> cov(std::complex<T> *vec, int comp, long int j) {
        // returns Sum_j C[i,j] vec[j]
        // Here, though, C0 is diagonal
        return C0s[comp][j]*vec[j+cumu_Ns[comp]];
    }

    //
    // GENERAL calculations that can run on constrained fields too
    //

    virtual std::complex<T> v1_dot_y(std::complex<T> *v1){
        // Calculate v1^{\dagger} y where y is the underlying realization
        T res_real(0), res_imag(0);


        for(int comp=0; comp<n_components; comp++) {
            #pragma omp parallel for reduction(+:res_real,res_imag)
            for(long int i=0; i<Ns[comp]; i++) {
                std::complex<T> res = conj(v1[i+cumu_Ns[comp]])*y0s[comp][i];
                // accumulate separately - OMP doesn't support complex number reduction :-(
                res_real+=std::real(res);
                res_imag+=std::imag(res);
            }
        }
        return std::complex<T>(res_real,res_imag);
    }

    virtual std::complex<T> v1_cov_v2(std::complex<T> *v1, std::complex<T> *v2){
        // Calculate v1^{\dagger} C v2
        T res_real(0), res_imag(0);
        for(int comp=0; comp<n_components; comp++) {
            #pragma omp parallel for reduction(+:res_real,res_imag)
            for(long int i=0; i<Ns[comp]; i++) {
                std::complex<T> res = conj(v1[i+cumu_Ns[comp]])*cov(v2,comp,i);
                // accumulate separately - OMP doesn't support complex number reduction :-(
                res_real+=std::real(res);
                res_imag+=std::imag(res);
            }
        }
        return std::complex<T>(res_real,res_imag);
    }

    virtual std::complex<T> v_cov_v(std::complex<T> *v){
        return v1_cov_v2(v,v);
    }


    void get_realization(std::complex<T> *r) {
        // Get the realization and store it in the specified array
        long int j=0;
        for(int comp=0; comp<n_components; comp++) {
            for(long int i=0; i<Ns[comp]; i++) {
                r[j]=y0s[comp][i];
                j++;
            }
        }
        assert(j==Ntot);
    }

    void get_mean(std::complex<T> *r) {
        // Get the mean of all realizations and store it in the specified array
        for(long int i=0; i<Ntot; i++) {
            if(i%1000==0) progress("get mean x", float(i)/Ntot);
            r[i]=0;
        }
        end_progress();
    }

    T get_delta_chi2() {
        return 0;
    }

};

template<typename T>
class MultiConstrainedField: public UnderlyingField<T>
{
private:

    std::vector<std::complex<T>* > alphas;
    std::vector<std::complex<T> > values;
    std::vector<std::complex<T> > existing_values;
public:
    UnderlyingField<T>* underlying;
    
    MultiConstrainedField(UnderlyingField<T>* underlying_, long int N) : UnderlyingField<T>(N),
    underlying(underlying_)
    {

    }

    void add_constraint(std::complex<T>* alpha, std::complex<T> value, std::complex<T> existing) {

        alphas.push_back(alpha);
        values.push_back(value);
        existing_values.push_back(existing);
    }

    void print_covariance() {
        int n=alphas.size();
        int done=0;
        std::vector<std::vector<std::complex<T>>> c_matrix(n,std::vector<std::complex<T>>(n,0));

        // Gram-Schmidt orthogonalization in-place
        for(int i=0; i<n; i++) {
            std::complex<T>* alpha_i=alphas[i];
            for(int j=0; j<=i; j++) {
                std::complex<T>* alpha_j=alphas[j];
                progress("calculating covariance", ((float)done*2)/(n*(1+n)));
                c_matrix[i][j]=underlying->v1_cov_v2(alpha_j,alpha_i);
                c_matrix[j][i]=c_matrix[i][j];
                done+=1;
            }
        }
        end_progress();

        std::cout << std::endl << "cov_matr = [";
        for(int i=0; i<n; i++) {
            std::cout << "[";
            for(int j=0; j<n; j++) {
                std::cout << std::real(c_matrix[i][j]);
                if(j<n-1) std::cout << ",";
            }
            std::cout << "]";
            if(i<n-1) std::cout << "," << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    void prepare() {

        int n=alphas.size();
        int done=0;

        std::cout << "v0=[";
        for(int i=0; i<n; i++) {
            std::cout << std::real(existing_values[i]);
            if(i<n-1)
                std::cout << ", ";
        }
        std::cout << "]"<<std::endl;

        std::cout << "v1=[";
        for(int i=0; i<n; i++) {
            std::cout << std::real(values[i]);
            if(i<n-1)
                std::cout << ", ";
            }
            std::cout << "]"<<std::endl;


        // Store transformation matrix, just for display purposes
        std::vector<std::vector<std::complex<T>>> t_matrix(n,std::vector<std::complex<T>>(n,0));

        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                if(i==j) t_matrix[i][j]=1.0; else t_matrix[i][j]=0.0;
            }
        }

        // Gram-Schmidt orthogonalization in-place
        for(int i=0; i<n; i++) {
            std::complex<T>* alpha_i=alphas[i];
            for(int j=0; j<i; j++) {
                std::complex<T>* alpha_j=alphas[j];
                progress("orthogonalizing constraints", ((float)done*2)/(n*(1+n)));
                std::complex<T> result = underlying->v1_cov_v2(alpha_j,alpha_i);

                #pragma omp parallel for
                for(long int k=0; k<this->Ntot; k++) {
                    alpha_i[k]-=result*alpha_j[k];
                }

                // update display matrix accordingly such that at each step
                // our current values array is equal to t_matrix . input_values
                for(int k=0; k<=j; k++)
                    t_matrix[i][k]-=result*t_matrix[j][k];

                // update constraining value
                values[i]-=result*values[j];
                done+=1; // one op for each of the orthognalizing constraints
            }

            // normalize
            std::complex<T> norm = sqrt(underlying->v_cov_v(alpha_i));
            #pragma omp parallel for
            for(long int k=0; k<this->Ntot; k++) {
                alpha_i[k]/=norm;
            }
            values[i]/=norm;

            for(int j=0;j<n; j++) {
                t_matrix[i][j]/=norm;
            }
            done+=1; // one op for the normalization

        }
        end_progress();

        // Now display t_matrix^{dagger} t_matrix, which is the matrix allowing one
        // to form chi^2: Delta chi^2 = d_1^dagger t_matrix^dagger t_matrix d_1 -
        // d_0^dagger t_matrix^dagger t_matrix d_0.



        existing_values.clear();
        for(int i=0; i<n; i++) {
            std::complex<T>* alpha_i=alphas[i];
            progress("calculating existing means",((float)i)/n);
            existing_values.push_back(underlying->v1_dot_y(alpha_i));
        }
        end_progress();

        std::cout << std::endl << "chi2_matr = [";
        for(int i=0; i<n; i++) {
            std::cout << "[";
            for(int j=0; j<n; j++) {
                std::complex<T> r=0;
                for(int k=0;k<n; k++) {
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

    T get_delta_chi2() {
        T rval=0.0, v0, v1;
        for(int i=0; i<alphas.size(); i++) {
            v0 = std::abs(existing_values[i]);
            v1 = std::abs(values[i]);
            rval+=v1*v1-v0*v0;
        }
        return rval/this->Ntot;
    }
    void get_realization(std::complex<T> *r) {
        int n=alphas.size();
        underlying->get_realization(r);

        for(int i=0; i<n; i++) {
            std::complex<T>* alpha_i=alphas[i];
            std::complex<T> dval_i = values[i]-existing_values[i];
            // std::cerr << "constraint " << i <<  values[i] << existing_values[i] << std::endl;
            progress("get constrained y", float(i)/n);
            #pragma omp parallel for
            for(long int j=0; j<this->Ntot; j++) {
                r[j]+=dval_i*underlying->cov(alpha_i,j);
            }
        }
        end_progress();
    }
};
