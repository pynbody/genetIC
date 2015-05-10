template<typename MyFloat>
std::complex<MyFloat> dot(std::complex<MyFloat> *a, std::complex<MyFloat> *b, const long int n) {
  std::complex<MyFloat> res(0);
  for(long int i=0; i<n; i++) res+=conj(a[i])*b[i];
  return res;
}

/*
template<typename MyFloat>
std::complex<MyFloat> dot(std::vector<std::complex<MyFloat>> &a,
                          std::vector<std::complex<MyFloat>> &b) {
  assert(a.size()==b.size());
  return dot(a.data(),b.data(),a.size());
}
*/


template<typename MyFloat>
std::complex<MyFloat> chi2(std::complex<MyFloat> *a, std::complex<MyFloat> *b, const long int n) {
    std::complex<MyFloat> res(0);
    for(long int i=1; i<n; i++) // go from 1 to avoid div by zero
    {
        res+=conj(a[i])*a[i]/b[i];
    }
    return res/((MyFloat)n);
}

/*
template<typename MyFloat>
std::complex<MyFloat> chi2(std::vector<std::complex<MyFloat>> &a,
                          std::vector<std::complex<MyFloat>> &b) {
  assert(a.size()==b.size());
  return chi2(a.data(),b.data(),a.size());
}
*/

template<typename MyFloat>
void mat_diag(const std::complex<MyFloat> *C, const std::complex<MyFloat> *alpha,
              const long int n, const std::complex<MyFloat> p,
              std::complex<MyFloat> *result)
{
  long int i;
  for(i=0;i<n;i++){ result[i]=C[i]*alpha[i]*p;}
}

template<typename MyFloat>
void mat_mat_diag(const std::complex<MyFloat> *z, const std::complex<MyFloat> *alpha, const long int n, const std::complex<MyFloat> p, std::complex<MyFloat> *result){
  long int i;
  std::complex<MyFloat> tf=0.;
  for(i=0;i<n;i++){
    tf+=z[i]*conj(alpha[i]);
  }
  for(i=0;i<n;i++){
    result[i]=alpha[i]*tf*p;
  }
}

template<typename MyFloat>
void mat_sum(const std::complex<MyFloat> *C1, const std::complex<MyFloat> *C2, const long int n, const std::complex<MyFloat> p, std::complex<MyFloat> *result){

  long int i;
  for(i=0;i<n;i++){
    result[i]=p*(C1[i]+C2[i]);
  }

}
