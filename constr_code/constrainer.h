
#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <Eigen/Dense>

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
  std::complex<T> *C0;
  std::complex<T> *y0;
  
protected:
  long int N;
  UnderlyingField(long int N): N(N) { }
  
public:
  UnderlyingField(std::complex<T> *C0,  std::complex<T> *y0, long int N) : C0(C0), 
									   y0(y0), N(N) {
  }

 
  //
  // SPECIFIC calculations for the underlying field
  //

  virtual std::complex<T> cov(std::complex<T> *vec, long int i) {
    // returns Sum_j C[i,j] vec[j]
    // Here, though, C0 is diagonal
    return C0[i]*vec[i];
  }

  //
  // GENERAL calculations that can run on constrained fields too
  //

   virtual std::complex<T> v1_dot_y(std::complex<T> *v1){
    // Calculate v1^{\dagger} y where y is the underlying realization
     T res_real(0), res_imag(0);
#pragma omp parallel for reduction(+:res_real,res_imag)
    for(long int i=0; i<N; i++) {
      std::complex<T> res = conj(v1[i])*y0[i];
      // accumulate separately - OMP doesn't support complex number reduction :-(
      res_real+=std::real(res);
      res_imag+=std::imag(res);
    }
    return std::complex<T>(res_real,res_imag);
  }

  virtual std::complex<T> v1_cov_v2(std::complex<T> *v1, std::complex<T> *v2){
    // Calculate v1^{\dagger} C v2
    T res_real(0), res_imag(0);
#pragma omp parallel for reduction(+:res_real,res_imag)
    for(long int i=0; i<N; i++) {
      std::complex<T> res = conj(v1[i])*cov(v2,i);
      // accumulate separately - OMP doesn't support complex number reduction :-(
      res_real+=std::real(res);
      res_imag+=std::imag(res);
    }
    return std::complex<T>(res_real,res_imag);
  }

  virtual std::complex<T> v_cov_v(std::complex<T> *v){
    return v1_cov_v2(v,v);
  }


  void get_realization(std::complex<T> *r) {
    // Get the realization and store it in the specified array
    for(long int i=0; i<N; i++) {
      r[i]=y0[i];
    }
  }

  void get_mean(std::complex<T> *r) {
    // Get the mean of all realizations and store it in the specified array
    for(long int i=0; i<N; i++) {
      if(i%1000==0) progress("get mean x", float(i)/N);
      r[i]=0;
    }
    end_progress();
  }

};

template<typename T>
class MultiConstrainedField: public UnderlyingField<T>
{
 private:
  UnderlyingField<T>* underlying;
  std::vector<std::complex<T>* > alphas;
  std::vector<std::complex<T> > values;
  std::vector<std::complex<T> > existing_values;
 public:
  MultiConstrainedField(UnderlyingField<T>* underlying_, long int N) : UnderlyingField<T>(N),
                                                                 underlying(underlying_)
   {
    
  }

  void add_constraint(std::complex<T>* alpha, std::complex<T> value) {
    
    alphas.push_back(alpha);
    values.push_back(value);
  }

  void prepare() {

    int n=alphas.size();
    int done=0;

    // Gram-Schmidt orthogonalization in-place
    for(int i=0; i<n; i++) {
      std::complex<T>* alpha_i=alphas[i];
      for(int j=0; j<i; j++) {
	std::complex<T>* alpha_j=alphas[j];
	progress("orthogonalizing constraints", ((float)done*2)/(n*(1+n)));
	std::complex<T> result = underlying->v1_cov_v2(alpha_j,alpha_i);
	#pragma omp parallel for
	for(long int k=0; k<this->N; k++) {
	  alpha_i[k]-=result*alpha_j[k];
	}
	values[i]-=result*values[j];
	done+=1; // one op for each of the orthognalizing constraints
      }

      // normalize
      std::complex<T> norm = sqrt(underlying->v_cov_v(alpha_i));
      #pragma omp parallel for
      for(long int k=0; k<this->N; k++) {
	alpha_i[k]/=norm;
      }
      values[i]/=norm;
      done+=1; // one op for the normalization

    }
    end_progress();

  
    
    existing_values.clear();
    for(int i=0; i<n; i++) {
      std::complex<T>* alpha_i=alphas[i];
      progress("calculating existing means",((float)i)/n);
      existing_values.push_back(underlying->v1_dot_y(alpha_i));
    }

    end_progress();
    
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
       for(long int j=0; j<this->N; j++) {
	 r[j]+=dval_i*underlying->cov(alpha_i,j);
       }
     }
    end_progress();
  }

   

};


