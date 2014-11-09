
#include <cmath>
#include <complex>
#include <map>


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
  std::cout << "                                                                                                       " << std::endl;
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

  virtual std::complex<T> mean(long int i) {
    // returns x0[i]
    return 0;
  }

  virtual std::complex<T> realization(long int i) {
    return y0[i];
  }

  virtual void destroy_stored_cov(std::complex<T> *v) {
    // only means something in the child class
  }

  //
  // GENERAL calculations that can run on constrained fields too
  //

  virtual std::complex<T> v_cov_v(std::complex<T> *v){
    // Calculate v^{\dagger} C v
    std::complex<T> res(0);
    for(long int i=0; i<N; i++) {
      res+=conj(v[i])*cov(v,i);
    }
    return res;
  }

  void get_realization(std::complex<T> *r) {
    // Get the realization and store it in the specified array
    for(long int i=0; i<N; i++) {
      if(i%1000==0) progress("get constrained y", float(i)/N);
      r[i]=realization(i);
    }
    end_progress();
  }

  void get_mean(std::complex<T> *r) {
    // Get the mean of all realizations and store it in the specified array
    for(long int i=0; i<N; i++) {
      if(i%1000==0) progress("get mean x", float(i)/N);
      r[i]=mean(i);
    }
    end_progress();
  }

};

template<typename T>
class ConstrainedField: public UnderlyingField<T>
{
private:
  UnderlyingField<T>* underlying;
  std::complex<T>* alpha;
  std::complex<T> alpha_dot_x0;
  std::complex<T> alpha_dot_d;
  std::complex<T> value;

  std::map<std::complex<T>*, std::complex<T> > cov_norm_values;

  std::complex<T> cov_norm(std::complex<T> *vec) {
    typename std::map<std::complex<T>*, std::complex<T> >::iterator search = cov_norm_values.find(vec);
    if(search==cov_norm_values.end()) {
      // not found - calcualte
      std::complex<T> value(0);
      for(long int i=0; i<this->N; i++)
	value+=std::conj(alpha[i])*underlying->cov(vec,i);
      cov_norm_values[vec] = value;
      return value;
    } else {
      return search->second;
    }

  }
public:
  
   ConstrainedField(UnderlyingField<T>* underlying_,
				std::complex<T>* alpha_,
				std::complex<T> value_,
		                long int N) : UnderlyingField<T>(N),
					      underlying(underlying_),
					      alpha(alpha_),
					      value(value_)
					                    
  {

    // normalize the constraint
    std::complex<T> norm(0);
    long int i;

    // normalize the constraint wrt C0
    for(i=0; i<N; i++)
      norm+=std::conj(alpha[i])*underlying->cov(alpha,i);

    norm = sqrt(norm);

    // renormalize alpha
    for(i=0; i<N; i++)
      alpha[i]/=norm;

    // we changed alpha so need to erase cached information on it
    underlying->destroy_stored_cov(alpha);
    
    // also renormalize the value
    value/=norm;
	
    // calculate alpha^dagger (y_0 - x_0) and alpha^dagger x_0    
    alpha_dot_d = 0;
    alpha_dot_x0 = 0;
    for(i=0; i<N; i++) {
      std::complex<T> x0_i = underlying->mean(i);
      alpha_dot_d+=std::conj(alpha[i])*(underlying->realization(i)-x0_i);
      alpha_dot_x0+=std::conj(alpha[i])*x0_i;
    }
  }

  virtual void destroy_stored_cov(std::complex<T> *v) {
    // forget any information we had about C0 v -- called if the contents
    // of v change
    cov_norm_values.erase(v);
    underlying->destroy_stored_cov(v);
  }

  virtual std::complex<T> cov(std::complex<T> *vec, long int i) {
    // Get C1_{ij} vec_j
    //
    // N.B. cov_norm term calculates alpha^{dagger} C0 vec, and caches it between
    // calls for efficiency
    std::complex<T> C0_vec_i = underlying->cov(vec,i);
    return C0_vec_i - underlying->cov(alpha,i)*cov_norm(vec);
  }

  virtual std::complex<T> mean(long int i) {
    // returns x1[i]
    return underlying->mean(i) + underlying->cov(alpha,i)*(value-alpha_dot_x0);
  }

  virtual std::complex<T> realization(long int i) {
    // returns y1[i]
    std::complex<T> d = underlying->realization(i) - underlying->mean(i);
    return d - alpha_dot_d*underlying->cov(alpha,i) + mean(i);
  }
  
};
