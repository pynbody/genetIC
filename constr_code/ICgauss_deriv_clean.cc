
#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "constrainer.h"

using namespace std;



#ifdef HAVE_HDF5
    #include "HDF_IO.hh"
#endif

#ifndef DOUBLEPRECISION     /* default is single-precision */
    typedef float  MyFloat;
    typedef float  MyDouble;
    #include <srfftw.h>

    #ifdef HAVE_HDF5
        hid_t hdf_float = H5Tcopy (H5T_NATIVE_FLOAT);
        hid_t hdf_double = H5Tcopy (H5T_NATIVE_FLOAT);
    #endif
    //#else
    //#if (DOUBLEPRECISION == 2)   /* mixed precision, do we want to implement this at all? */
    //typedef float   MyFloat;
    //typedef double  MyDouble;
    //hid_t hdf_float = H5Tcopy (H5T_NATIVE_FLOAT);
    //hid_t hdf_double = H5Tcopy (H5T_NATIVE_DOUBLE);
#else                        /* everything double-precision */
    typedef double  MyFloat;
    typedef double  MyDouble;


    #ifdef FFTW_TYPE_PREFIX
    #include <drfftw.h>
    #else
    #include <rfftw.h>
    #endif
#endif




typedef struct _ic_state {
    int *part_arr;
    int n_part_arr;
    MyFloat boxlen;
    int res;
    long nPartTotal;
    MyFloat dx; // = boxlen/res
    MyFloat x0,y0,z0;
    float a, om;
    Grid<MyFloat> *pGrid;
    std::complex<MyFloat> *pField_k;
    std::complex<MyFloat> *pField_x;
} ic_state;


MyFloat get_wrapped_delta(MyFloat x0, MyFloat x1, const ic_state *pState) {
    MyFloat boxlen = pState->boxlen;
    MyFloat result = x0-x1;
    if(result>boxlen/2) {
        result-=boxlen;
    }
    if(result<-boxlen/2) {
        result+=boxlen;
    }
    return result;
}


void getCentre(ic_state *pIcs) {

    pIcs->x0 = 0; pIcs->y0 = 0; pIcs->z0 =0;

    MyFloat xa=pIcs->pGrid->cells[pIcs->part_arr[0]].coords[0]*pIcs->dx;
    MyFloat ya=pIcs->pGrid->cells[pIcs->part_arr[0]].coords[1]*pIcs->dx;
    MyFloat za=pIcs->pGrid->cells[pIcs->part_arr[0]].coords[2]*pIcs->dx;

    for(long i=0;i<pIcs->n_part_arr;i++) {
        pIcs->x0+=get_wrapped_delta(pIcs->pGrid->cells[pIcs->part_arr[i]].coords[0]*pIcs->dx,xa,pIcs);
        pIcs->y0+=get_wrapped_delta(pIcs->pGrid->cells[pIcs->part_arr[i]].coords[1]*pIcs->dx,ya,pIcs);
        pIcs->z0+=get_wrapped_delta(pIcs->pGrid->cells[pIcs->part_arr[i]].coords[2]*pIcs->dx,za,pIcs);
    }
    pIcs->x0/=pIcs->n_part_arr;
    pIcs->y0/=pIcs->n_part_arr;
    pIcs->z0/=pIcs->n_part_arr;
    pIcs->x0+=xa;
    pIcs->y0+=ya;
    pIcs->z0+=za;
}

extern void ensureRealDelta(ic_state *pIcs);

void centreDenmax(ic_state *pIcs) {

    float den_max=-1000;
    long index_max=0;
    ensureRealDelta(pIcs);

    for(long i=0;i<pIcs->n_part_arr;i++) {
        if(std::real(pIcs->pField_x[pIcs->part_arr[i]])>den_max) {
            index_max=pIcs->part_arr[i];
            den_max = std::real(pIcs->pField_x[pIcs->part_arr[i]]);
            pIcs->x0 = pIcs->pGrid->cells[pIcs->part_arr[i]].coords[0]*pIcs->dx;
            pIcs->y0 = pIcs->pGrid->cells[pIcs->part_arr[i]].coords[1]*pIcs->dx;
            pIcs->z0 = pIcs->pGrid->cells[pIcs->part_arr[i]].coords[2]*pIcs->dx;
        }
    }

    cerr << "Denmax = " << den_max <<", index=" << index_max << " coords=" << pIcs->x0 << " " << pIcs->y0 << " " << pIcs->z0 << endl;

}

void centreOn(ic_state *pIcs, long id) {
    pIcs->x0 = pIcs->pGrid->cells[id].coords[0]*pIcs->dx;
    pIcs->y0 = pIcs->pGrid->cells[id].coords[1]*pIcs->dx;
    pIcs->z0 = pIcs->pGrid->cells[id].coords[2]*pIcs->dx;
}

void selectSphere(ic_state *pIcs, float radius) {
    float r2 = radius*radius;
    float delta_x, delta_y, delta_z, r2_i;
    int n=0;
    for(long i=0;i<pIcs->nPartTotal;i++) {
        delta_x = get_wrapped_delta(pIcs->pGrid->cells[i].coords[0]*pIcs->dx,pIcs->x0,pIcs);
        delta_y = get_wrapped_delta(pIcs->pGrid->cells[i].coords[1]*pIcs->dx,pIcs->y0,pIcs);
        delta_z = get_wrapped_delta(pIcs->pGrid->cells[i].coords[2]*pIcs->dx,pIcs->z0,pIcs);
        r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
        if(r2_i<r2)
            n++;
    }
    cerr << "Selecting " << n << " particles..." << endl;
    if(pIcs->part_arr!=NULL)
        free(pIcs->part_arr);

    pIcs->part_arr = (int*)calloc(n,sizeof(int));
    n=0;
    for(long i=0;i<pIcs->nPartTotal;i++) {
        delta_x = get_wrapped_delta(pIcs->pGrid->cells[i].coords[0]*pIcs->dx,pIcs->x0,pIcs);
        delta_y = get_wrapped_delta(pIcs->pGrid->cells[i].coords[1]*pIcs->dx,pIcs->y0,pIcs);
        delta_z = get_wrapped_delta(pIcs->pGrid->cells[i].coords[2]*pIcs->dx,pIcs->z0,pIcs);
        r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
        if(r2_i<r2)
            pIcs->part_arr[n++]=i;
    }
    pIcs->n_part_arr = n;


}
void AllocAndGetBuffer_int(const char* IDfile, ic_state *pIcs, bool append=false) {
    // count lines, allocate space, then read

    FILE *f;
    int i=0;
    int r=0;
    int c;

    f = fopen(IDfile, "r");
    if (f == NULL) {printf("\n Input file %s not found!\n", IDfile); exit(-1);}

    while ( (c=fgetc(f)) != EOF ) {
        if ( c == '\n' )
                i++;
    }

    fclose(f);

    cerr << "File " <<IDfile << " has " << i << " lines" << endl;

    if(append && pIcs->part_arr==NULL) {
        cerr << "Can't append - no particles were loaded" << endl;
        append=false;
    }
    int final_size = i;

    if(append) final_size+=pIcs->n_part_arr;
    int *arr = (int*)calloc(final_size,sizeof(int));

    GetBuffer_int(arr, IDfile, i);

    // copy old particles into new buffer
    if (append>0)
      cerr << "Also keeping " << pIcs->n_part_arr << " existing particles" << endl;
    for(int j=0; j<pIcs->n_part_arr; j++)
      arr[i+j] = pIcs->part_arr[j];

    if(pIcs->part_arr!=NULL)
      free(pIcs->part_arr);

    pIcs->part_arr = arr;
    pIcs->n_part_arr = final_size;

    getCentre(pIcs);

}

string make_base( string basename, int res, MyFloat box, MyFloat zin){
  ostringstream result;
#ifndef DOUBLEPRECISION
  result << basename<<"IC_iter_sing_z"<<zin<<"_"<<res<<"_L" << box;
#else
  result << basename<<"IC_iter_doub_z"<<zin<<"_"<<res<<"_L"<< box;
#endif

  return result.str();

}

MyFloat sig(MyFloat R, double *kcamb, double *Tcamb, MyFloat ns, MyFloat L, int res, int quoppa) {

  MyFloat s=0.,k,t;

  MyFloat amp=9./2./M_PI/M_PI;
  MyFloat kmax=kcamb[quoppa-1];
  MyFloat kmin=kcamb[0];

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, quoppa);
  gsl_spline_init (spline, kcamb, Tcamb, quoppa);

  MyFloat dk=(kmax-kmin)/10000.;
    for(k=kmin; k<kmax;k+=dk){

	t=gsl_spline_eval (spline, k, acc);

	s+= powf(k,ns+2.)*((sin(k*R)-k*R*cos(k*R))/((k*R)*(k*R)*(k*R)) ) *( (sin(k*R)-k*R*cos(k*R))/((k*R)*(k*R)*(k*R)) )*t*t;

    }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  s=sqrt(s*amp*dk);
  return s;

}

MyFloat D(MyFloat a, MyFloat Om, MyFloat Ol){

   MyFloat Hsq=Om/powf(a,3.0)+(1.-Om-Ol)/a/a+Ol;
   MyFloat d=2.5*a*Om/powf(a,3.0)/Hsq/( powf(Om/Hsq/a/a/a,4./7.) - Ol/Hsq + (1.+0.5*Om/powf(a,3.0)/Hsq)*(1.+1./70.*Ol/Hsq) );

   //simplify this...?

 return d;
}


complex<MyFloat> *fft_r(complex<MyFloat> *fto, complex<MyFloat> *ftin, const int res, const int dir){ //works when fto and ftin and output are the same, but overwrites fto for the backwards transform so we have another temporary array there

    long i;

  if(dir==1){
	std::complex<MyFloat> *fti = (std::complex<MyFloat>*)calloc(res*res*res,sizeof(std::complex<MyFloat>));
	for(i=0;i<res*res*res;i++){fti[i]=ftin[i]/sqrt((MyFloat)(res*res*res));}
	//for(i=0;i<res*res*res;i++){fti[i]=ftin[i]/(MyFloat)(res*res*res);}
	fftwnd_plan pf = fftw3d_create_plan( res, res, res, FFTW_FORWARD, FFTW_ESTIMATE );
	fftwnd_one(pf, reinterpret_cast<fftw_complex*>(&fti[0]), reinterpret_cast<fftw_complex*>(&fto[0]));
	fftwnd_destroy_plan(pf);
	free(fti);
	    }

  else if(dir==-1){
	std::complex<MyFloat> *fti = (std::complex<MyFloat>*)calloc(res*res*res,sizeof(std::complex<MyFloat>));
	//memcpy(fti,ftin,res*res*res*sizeof(fftw_complex));
	for(i=0;i<res*res*res;i++){fti[i]=ftin[i]/sqrt((MyFloat)(res*res*res));}
	fftwnd_plan pb;
 	pb = fftw3d_create_plan( res,res,res, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftwnd_one(pb, reinterpret_cast<fftw_complex*>(&fti[0]), reinterpret_cast<fftw_complex*>(&fto[0]));
	fftwnd_destroy_plan(pb);
	free(fti);
	  }

  else {std::cerr<<"wrong parameter for direction in fft_r"<<std::endl;}

  return fto;
}

complex<MyFloat> dot(complex<MyFloat> *a, complex<MyFloat> *b, const long int n) {
  complex<MyFloat> res(0);
  for(long int i=0; i<n; i++) res+=conj(a[i])*b[i];
  return res;
}

complex<MyFloat> chi2(complex<MyFloat> *a, complex<MyFloat> *b, const long int n) {
    complex<MyFloat> res(0);
    for(long int i=1; i<n; i++) // go from 1 to avoid div by zero
    {
        res+=conj(a[i])*a[i]/b[i];
    }
    return res/((MyFloat)n);
}

void mat_diag(const complex<MyFloat> *C, const complex<MyFloat> *alpha, const long int n, const complex<MyFloat> p, complex<MyFloat> *result){

  long int i;
  for(i=0;i<n;i++){ result[i]=C[i]*alpha[i]*p;}

}

void mat_mat_diag(const complex<MyFloat> *z, const complex<MyFloat> *alpha, const long int n, const complex<MyFloat> p, complex<MyFloat> *result){

  long int i;
  complex<MyFloat> tf=0.;
  for(i=0;i<n;i++){
    tf+=z[i]*conj(alpha[i]);
  }
  for(i=0;i<n;i++){
    result[i]=alpha[i]*tf*p;
  }


}

void mat_sum(const complex<MyFloat> *C1, const complex<MyFloat> *C2, const long int n, const complex<MyFloat> p, complex<MyFloat> *result){

  long int i;
  for(i=0;i<n;i++){
    result[i]=p*(C1[i]+C2[i]);
  }

}


complex<MyFloat>* calc_A_new(complex<MyFloat>* z, const complex<MyFloat> *alpha, const long int n, const complex<MyFloat> *C0, complex<MyFloat> *temp, complex<MyFloat> *temp2, complex<MyFloat> *temp3){
  // Added by AP - non-iterative version... sorry! :-/
  mat_mat_diag(z,alpha,n, complex<MyFloat>(-1.0,0.), temp);
  return temp;
}

complex<MyFloat>* calc_A(int iter, complex<MyFloat>* z, const complex<MyFloat> *alpha, const long int n, const complex<MyFloat> *C0, complex<MyFloat> *temp, complex<MyFloat> *temp2, complex<MyFloat> *temp3){

  if(iter==0){
    mat_mat_diag(z,alpha,n, complex<MyFloat>(-0.5,0.), temp);
    return temp;
  }
  else{

  temp3=calc_A(iter-1,z, alpha, n, C0, temp, temp2, temp3);

  complex<MyFloat> *md=(complex<MyFloat>*)calloc(n,sizeof(complex<MyFloat>));

  mat_diag(C0,temp3, n, complex<MyFloat>(1.,0.), md);

  temp2=calc_A(iter-1, md, alpha, n, C0, temp, temp3, temp2);

  mat_mat_diag(z,alpha,n, complex<MyFloat>(1.,0.), md);
  mat_sum(temp2, md, n, complex<MyFloat>(-0.5,0.), temp3);

  free(md);

  return temp3;

  }
}


void ret_exp(long res, MyFloat R, MyFloat mr, MyFloat L, std::complex<MyFloat> *ret){//generates a Gaussian constraint vector

  long r1, r2, r3, i;//, rr1, rr2, rr3;
  MyFloat r, w, n_r=0.;

  MyFloat rw=1./MyFloat(res)*L;

        for(r1=0; r1<res;r1++){
  	  for(r2=0;r2<res;r2++){
      	    for(r3=0;r3<res;r3++){

		i = (r1*res+r2)*(res)+r3;
		r=MyFloat( ( (r1-mr)*(r1-mr) + (r2-mr)*(r2-mr) + (r3-mr)*(r3-mr) ) );
		w=(MyFloat)exp(-(r/R/R*rw*rw/2.));

		ret[i].real(w);
		n_r+=w;

	    }
	  }
	}

     for(i=0; i<res*res*res; i++){ ret[i]/=n_r;}

}


//faster than sorting, even though we interpolate more often
void brute_interpol_new(int res, double *kcamb, double *Tcamb, int quoppa, MyFloat kw, MyFloat ns, MyFloat norm_amp, complex<MyFloat> *ft, complex<MyFloat> *ftsc, complex<MyFloat> *P){

  complex<MyFloat> norm_iter (0.,0.);

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, quoppa);
  gsl_spline_init (spline, kcamb, Tcamb, quoppa);

  MyFloat kk=0.;
  long k1,k2,k3;
  long ii,jj,ll;
  long idk;

  for(k1=-res/2; k1<res/2;k1++){
      if(k1<0){ii=k1+res;}else{ii=k1;}
     for(k2=-res/2;k2<res/2;k2++){
       if(k2<0){jj=k2+res;}else{jj=k2;}
       for(k3=-res/2;k3<res/2;k3++){
	 if(k3<0){ll=k3+res;}else{ll=k3;}


	 idk=(ii*(res)+jj)*(res)+ll;
	 if(idk==0){continue;}
	 kk=sqrt(k1*k1+k2*k2+k3*k3)*kw;

	 P[idk]=gsl_spline_eval (spline, kk, acc);
	 P[idk]*=(P[idk]* powf(kk,ns) *norm_amp );

	 ftsc[idk]=sqrt(P[idk])*ft[idk];

       }
    }
  }


  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

}


void powsp(int n, complex<MyFloat>* ft, const char *out, MyFloat Boxlength){

 int res=n;
 int nBins=100;
 MyFloat *inBin = new MyFloat[nBins];
 MyFloat *Gx = new MyFloat[nBins];
 MyFloat *kbin = new MyFloat[nBins];
 MyFloat kmax = M_PI/Boxlength*(MyFloat)res, kmin = 2.0f*M_PI/(MyFloat)Boxlength, dklog = log10(kmax/kmin)/nBins, kw = 2.0f*M_PI/(MyFloat)Boxlength;

      int ix,iy,iz,idx,idx2;
      MyFloat kfft;

      for( ix=0; ix<nBins; ix++ ){
	inBin[ix] = 0;
	Gx[ix] = 0.0f;
	kbin[ix] = 0.0f;
      }


   for(ix=0; ix<res;ix++)
    for(iy=0;iy<res;iy++)
      for(iz=0;iz<res;iz++){
        idx = (ix*res+iy)*(res)+iz;

        // determine mode modulus

	MyFloat vabs = ft[idx].real()*ft[idx].real()+ft[idx].imag()*ft[idx].imag();

	int iix, iiy, iiz;

        if( ix>res/2 ) iix = ix - res; else iix = ix;
        if( iy>res/2 ) iiy = iy - res; else iiy = iy;
        if( iz>res/2 ) iiz = iz - res; else iiz = iz;

        kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);
       MyFloat k = kfft*kw;

        // correct for aliasing, formula from Jing (2005), ApJ 620, 559
        // assume isotropic aliasing (approx. true for k<kmax=knyquist)
        // this formula is for CIC interpolation scheme, which we use <- only needed for Powerspectrum

       MyFloat JingCorr = (1.0f-2.0f/3.0f*sin(M_PI*k/kmax/2.0f)*sin(M_PI*k/kmax/2.0f));
	vabs /= JingCorr;

        //.. logarithmic spacing in k
	      idx2 = (int)(( 1.0f/dklog * log10(k/kmin) )); //cout<<idx2<<endl;

        if(k>=kmin&&k<kmax){


            Gx[idx2] += vabs;
            kbin[idx2] += k;
            inBin[idx2]++;

        }else{continue;}

      }


  //... convert to physical units ...
  std::ofstream ofs(out);

  //definition of powerspectrum brings (2pi)^-3, FT+conversion to physical units brings sqrt(Box^3/N^6) per delta1, where ps22~d2*d2~d1*d1*d1*d1 -> (Box^3/N^6)^2

	MyFloat psnorm=powf(Boxlength/(2.0*M_PI),3.0);

  for(ix=0; ix<nBins; ix++){

	if(inBin[ix]>0){


    ofs << std::setw(16) << pow (10., log10 (kmin) + dklog * (ix + 0.5) )
	<< std::setw(16) <<  kbin[ix]/inBin[ix]
        << std::setw(16) << (MyFloat)(Gx[ix]/inBin[ix])*psnorm
        << std::setw(16) << inBin[ix]
        << std::endl;

	}
  }
  ofs.close();

  free(inBin);
  free(kbin);
  free(Gx);

}

void powsp_noJing(int n, complex<MyFloat>* ft, const char *out, MyFloat Boxlength){

 int res=n;
 int nBins=100;
 MyFloat *inBin = new MyFloat[nBins];
 MyFloat *Gx = new MyFloat[nBins];
 MyFloat *kbin = new MyFloat[nBins];
 MyFloat kmax = M_PI/Boxlength*(MyFloat)res, kmin = 2.0f*M_PI/(MyFloat)Boxlength, dklog = log10(kmax/kmin)/nBins, kw = 2.0f*M_PI/(MyFloat)Boxlength;

      int ix,iy,iz,idx,idx2;
      MyFloat kfft;

      for( ix=0; ix<nBins; ix++ ){
	inBin[ix] = 0;
	Gx[ix] = 0.0f;
	kbin[ix] = 0.0f;
      }


   for(ix=0; ix<res;ix++)
    for(iy=0;iy<res;iy++)
      for(iz=0;iz<res;iz++){
        idx = (ix*res+iy)*(res)+iz;

        // determine mode modulus

	MyFloat vabs = ft[idx].real()*ft[idx].real()+ft[idx].imag()*ft[idx].imag();

	int iix, iiy, iiz;

        if( ix>res/2 ) iix = ix - res; else iix = ix;
        if( iy>res/2 ) iiy = iy - res; else iiy = iy;
        if( iz>res/2 ) iiz = iz - res; else iiz = iz;

        kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);
       MyFloat k = kfft*kw;

        //.. logarithmic spacing in k
	      idx2 = (int)(( 1.0f/dklog * log10(k/kmin) ));

        if(k>=kmin&&k<kmax){

            Gx[idx2] += vabs/(MyFloat)(res*res*res); //because FFT is now normalised with 1/sqrt(Ntot)
            kbin[idx2] += k;
            inBin[idx2]++;

        }else{continue;}

      }


  //... convert to physical units ...
  std::ofstream ofs(out);


	MyFloat psnorm=powf(Boxlength/(2.0*M_PI),3.0);

  for(ix=0; ix<nBins; ix++){

	if(inBin[ix]>0){


    ofs << std::setw(16) << pow (10., log10 (kmin) + dklog * (ix + 0.5) )
	<< std::setw(16) <<  kbin[ix]/inBin[ix]
        << std::setw(16) << (MyFloat)(Gx[ix]/inBin[ix])*psnorm
        << std::setw(16) << inBin[ix]
        << std::endl;

	}
  }
  ofs.close();

  free(inBin);
  free(kbin);
  free(Gx);

}

complex<MyFloat> *poiss(complex<MyFloat>* out, complex<MyFloat> *in, int res, MyFloat Boxlength, MyFloat a, MyFloat Om){

  long i;
  MyFloat prefac=3./2.*Om/a*100.*100./(3.*100000.)/(3.*100000.); // =3/2 Om0/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)
  MyFloat  kw = 2.0f*M_PI/(MyFloat)Boxlength,k;

  int k1,k2,k3,kk1,kk2,kk3;

   for(k1=0; k1<res;k1++){
       for(k2=0;k2<res;k2++){
          for(k3=0;k3<res;k3++){
        	i = (k1*res+k2)*(res)+k3;

		if( k1>res/2 ) kk1 = k1-res; else kk1 = k1;
        	if( k2>res/2 ) kk2 = k2-res; else kk2 = k2;
 		if( k3>res/2 ) kk3 = k3-res; else kk3 = k3;

		k=(MyFloat)(kk1*kk1+kk2*kk2+kk3*kk3)*kw*kw;

		out[i]=-in[i]*prefac/k;
	  }
       }
   }

   out[0]=complex<MyFloat>(0.,0.);

   return out;

}

complex<MyFloat> *rev_poiss(complex<MyFloat> *out, complex<MyFloat> *in, int res, MyFloat Boxlength, MyFloat a, MyFloat Om){

  long i;

  MyFloat prefac=3./2.*Om/a*100.*100./(3.*100000.)/(3.*100000.); // 3/2 Om/a * (H0/h)^2 (h/Mpc)^2 / c^2 (km/s)

  MyFloat  kw = 2.0f*M_PI/(MyFloat)(Boxlength),k;

  int k1,k2,k3,kk1,kk2,kk3;

   for(k1=0; k1<res;k1++){
      for(k2=0;k2<res;k2++){
         for(k3=0;k3<res;k3++){
	    i = (k1*res+k2)*(res)+k3;

		if( k1>res/2 ) kk1 = k1-res; else kk1 = k1;
        	if( k2>res/2 ) kk2 = k2-res; else kk2 = k2;
 		if( k3>res/2 ) kk3 = k3-res; else kk3 = k3;

		k=(MyFloat)(kk1*kk1+kk2*kk2+kk3*kk3)*kw*kw;

		out[i]=-k*in[i]/prefac;
	 }
      }
   }

    out[0]=complex<MyFloat> (0.,0.);

   return out;

}

void pbc(long n, std::complex<MyFloat> *delta_in, std::complex<MyFloat> *delta_out, long r1, long r2, long r3, long *index_shift){

  if(abs(r1) > n ||  abs(r2) > n || abs(r3) > n ){std::cerr << "Wrong shift parameters: " << r1 << " "<< r2 << " "<<r3 <<std::endl; exit(1);}

  long ix,iy,iz,id1,id2,jx,jy,jz,jjx,jjy,jjz;
  for(ix=0;ix<n;ix++){
    for(iy=0;iy<n;iy++){
      for(iz=0;iz<n;iz++){

	       id1=(ix*n+iy)*n+iz;

	       jx=ix-r1;
	       jy=iy-r2;
	       jz=iz-r3;

	       if( jx>n-1 ) jjx = jx - n; else if(jx <0) jjx = jx + n; else jjx = jx;
	       if( jy>n-1 ) jjy = jy - n; else if(jy <0) jjy = jy + n; else jjy = jy;
	       if( jz>n-1 ) jjz = jz - n; else if(jz <0) jjz = jz + n; else jjz = jz;

	       id2=(jjx*n+jjy)*n+jjz;

	       delta_out[id2]=delta_in[id1];
         index_shift[id1]=id2;

      }
    }
  }


}

void GridParticles( int nx, int ny, int nz, MyFloat *Pos1,MyFloat *Pos2,MyFloat *Pos3, MyFloat Boxlength, unsigned nPart, MyFloat m,  fftw_real *data , MyFloat subbox, int *ourbox){

  unsigned ix, iy, iz, i;

  MyFloat wpar =  (MyFloat)(nx*ny*nz)/(MyFloat)nPart;

  MyFloat x,y,z,dx,dy,dz,tx,ty,tz,tyw,dyw;
  unsigned ix1,iy1,iz1;
  int icount = 0;
  unsigned nPartThis = nPart;

  for (i = 0; i < 3; i++)
      if (ourbox[i] > subbox || ourbox[i] < 1){cerr<< "Error with ourbox."; exit(-1);}
	  //error (-1,0,"Error with ourbox.");

  for( unsigned i=0; i<nPartThis; ++i ){
    x = Pos1[i];
    y = Pos2[i];
    z = Pos3[i];



    x -= (ourbox[0]-1.) * (MyFloat)nx;
    y -= (ourbox[1]-1.) * (MyFloat)ny;
    z -= (ourbox[2]-1.) * (MyFloat)nz;

      if (x < nx && y < ny && z < nz && x > 0. && y > 0. && z > 0.){
	  ++icount;
      }
  }

  wpar =  (MyFloat)(nx*ny*nz)/(MyFloat)nPart;


  for( unsigned i=0; i<nPartThis; ++i ){

    x = Pos1[i];
    y = Pos2[i];
    z = Pos3[i];


    x -= (ourbox[0]-1.) * (MyFloat)nx;
    y -= (ourbox[1]-1.) * (MyFloat)ny;
    z -= (ourbox[2]-1.) * (MyFloat)nz;


    if (x < nx && y < ny && z < nz && x > 0. && y > 0. && z > 0.){


    ix = (unsigned)x;
    iy = (unsigned)y;
    iz = (unsigned)z;

    dx = (x-((MyFloat)ix));
    dy = (y-((MyFloat)iy));
    dz = (z-((MyFloat)iz));

    ix %= nx;
    iy %= ny;
    iz %= nz;

    tx = 1.0f-dx;
    ty = 1.0f-dy;
    tz = 1.0f-dz;

    tyw = ty*wpar;
    dyw = dy*wpar;

    ix1 = (ix+1)%nx;
    iy1 = (iy+1)%ny;
    iz1 = (iz+1)%nz;

    data[(ix*ny + iy) * nz + iz]   += tz*tx*tyw;
    data[(ix1*ny + iy) * nz + iz]  += tz*dx*tyw;
    data[(ix*ny + iy1) * nz + iz]  += tz*tx*dyw;
    data[(ix1*ny + iy1) * nz + iz] += tz*dx*dyw;

    data[(ix*ny + iy) * nz + iz1]   += dz*tx*tyw;
    data[(ix1*ny + iy) * nz + iz1]  += dz*dx*tyw;
    data[(ix*ny + iy1) * nz + iz1]  += dz*tx*dyw;
    data[(ix1*ny + iy1) * nz + iz1] += dz*dx*dyw;
    }

  }
}

void old_check_pbc_coords(long *coords, int size){

  while(coords[0]> size){coords[0]-=size;}
  while(coords[1]> size){coords[1]-=size;}
  while(coords[2]> size){coords[2]-=size;}

  while(coords[0]<0){coords[0]+=size;}
  while(coords[1]<0){coords[1]+=size;}
  while(coords[2]<0){coords[2]+=size;}


}



void cen_deriv4_alpha(long index, int direc, complex<MyFloat> *alpha, const ic_state *pState)

{//4th order central difference

  MyFloat x0, y0, z0;//, zm2, zm1, zp2, zp1;
  x0=get_wrapped_delta(pState->dx*pState->pGrid->cells[index].coords[0],pState->x0,pState);
  y0=get_wrapped_delta(pState->dx*pState->pGrid->cells[index].coords[1],pState->y0,pState);
  z0=get_wrapped_delta(pState->dx*pState->pGrid->cells[index].coords[2],pState->z0,pState);

  complex<MyFloat> ang=0.;

  //do the epsilon
  MyFloat c[3]={0,0,0};
  if(direc==0){c[2]=y0; c[1]=-z0;}
  else if(direc==1){c[0]=z0; c[2]=-x0;}
  else if(direc==2){c[1]=x0; c[0]=-y0;}
  else if(direc==3){
      MyFloat r0 = std::sqrt((x0*x0)+(y0*y0)+(z0*z0));
      if(r0!=0) {
          c[0]=x0/r0;
          c[1]=y0/r0;
          c[2]=z0/r0;
      }
  } // radial velocity

  else{cerr<< "Wrong value for parameter 'direc' in function 'cen_deriv4_alpha'."<< endl; exit(1);}


  for(int di=0; di<3; di++) {
      long ind_p1, ind_m1, ind_p2, ind_m2;
      //first step in rho direction
      int step1[3]={0,0,0};
      int neg_step1[3]={0,0,0};
      step1[di]=1;
      neg_step1[di]=-1;

      ind_m1=pState->pGrid->find_next_ind(index, neg_step1);
      ind_p1=pState->pGrid->find_next_ind(index, step1);


      ind_m2=pState->pGrid->find_next_ind(ind_m1, neg_step1);
      ind_p2=pState->pGrid->find_next_ind(ind_p1, step1);

      MyFloat a=-1./12./pState->dx, b=2./3./pState->dx;  //the signs here so that L ~ - Nabla Phi

      alpha[ind_m2]+=(c[di]*a);
      alpha[ind_m1]+=(c[di]*b);
      alpha[ind_p1]+=(-c[di]*b);
      alpha[ind_p2]+=(-c[di]*a);
   }

}

/*
complex<MyFloat> cen_deriv2_alpha(Grid *in, int index, complex<MyFloat> *alpha, MyFloat dx, int direc, complex<MyFloat> *phi){//second order central difference

  MyFloat x0, y0, z0;
  x0=in->cells[index].coords[0];
  y0=in->cells[index].coords[1];
  z0=in->cells[index].coords[2];

  complex<MyFloat> ang=0.;

  //do the epsilon
  MyFloat c1,c2;
  int d1,d2;
  if(direc==0){d2=1; d1=2; c1=y0; c2=z0;}
  else if(direc==1){d2=2; d1=0, c1=z0; c2=x0;}
  else if(direc==2){d2=0; d1=1; c1=x0; c2=y0;}
  else{cerr<< "Wrong value for parameter 'direc' in function 'cen_deriv2_alpha'."<< endl; exit(1);}

  int ind_p1, ind_m1;
  //first step in rho direction
    int step1[3]={0,0,0};
    step1[d1]=1;
    int neg_step1[3]={0,0,0};
    neg_step1[d1]=-1;

    ind_m1=in->find_next_ind(index, neg_step1);
    ind_p1=in->find_next_ind(index, step1);

    MyFloat a=-0.5/dx;

    alpha[ind_m1]+=(c1*a);
    alpha[ind_p1]+=(-c1*a);

    ang+=((c1*a*phi[ind_m1])+(-c1*a*phi[ind_p1]));

  //second step in other rho direction
    int step2[3]={0,0,0};
    step2[d2]=1;
    int neg_step2[3]={0,0,0};
    neg_step2[d2]=-1;

    ind_m1=in->find_next_ind(index, neg_step2);
    ind_p1=in->find_next_ind(index, step2);

    alpha[ind_m1]+=(-c2*a);
    alpha[ind_p1]+=(c2*a);

    ang+=((-c2*a*phi[ind_m1])+(c2*a*phi[ind_p1]));

    return ang;

}

*/

void reorderBuffer(ic_state *pIcs) {

    cout << "Reordering buffer radially..." << endl;
    cout << " [taking centre = " << pIcs->x0 << " " << pIcs->y0 << " " <<pIcs->z0 << "]" << endl;

    int *part_arr = pIcs->part_arr;
    int n_part_arr = pIcs->n_part_arr;

    MyFloat r2[n_part_arr];
    MyFloat x0=pIcs->x0,y0=pIcs->y0,z0=pIcs->z0;
    MyFloat delta_x, delta_y, delta_z;
    std::vector<size_t> index(n_part_arr);

    for(int i=0;i<n_part_arr;i++) {
        delta_x = get_wrapped_delta(pIcs->pGrid->cells[part_arr[i]].coords[0]*pIcs->dx,x0,pIcs);
        delta_y = get_wrapped_delta(pIcs->pGrid->cells[part_arr[i]].coords[1]*pIcs->dx,y0,pIcs);
        delta_z = get_wrapped_delta(pIcs->pGrid->cells[part_arr[i]].coords[2]*pIcs->dx,z0,pIcs);
        r2[i] = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
        index[i]=i;
    }

    // Now sort the index array
    std::sort(index.begin(),index.end(),
         [&r2](size_t i1, size_t i2) { return r2[i1]<r2[i2]; } );

    // Turn the index array into something pointing to the particles
    for(int i=0; i<n_part_arr; i++) {
        index[i] = part_arr[index[i]];
    }

    // Copy back into the particle array
    for(int i=0; i<n_part_arr; i++) {
        part_arr[i] = index[i];
    }

}

complex<MyFloat> *calcConstraintVector(istream &inf, const ic_state *pState) {

    char name[100];
    inf >> name;
    cout << "Getting constraint vector '" << name << "' for " << pState->n_part_arr << " particles.";

    complex<MyFloat> *rval=(complex<MyFloat>*)calloc(pState->nPartTotal,sizeof(complex<MyFloat>));
    complex<MyFloat> *rval_k=(complex<MyFloat>*)calloc(pState->nPartTotal,sizeof(complex<MyFloat>));

    if(strcasecmp(name,"overdensity")==0) {
        MyFloat w = 1.0/pState->n_part_arr;
        for(long i=0;i<pState->n_part_arr;i++) {
            rval[pState->part_arr[i]]=w;
        }

        fft_r(rval_k, rval, pState->res, 1);
    }
    else if(strcasecmp(name,"phi")==0) {
        MyFloat w = 1.0/pState->n_part_arr;
        for(long i=0;i<pState->n_part_arr;i++) {
            rval[pState->part_arr[i]]=w;
        }
        complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(pState->nPartTotal,sizeof(complex<MyFloat>));
        fft_r(rval_kX, rval, pState->res, 1);
        poiss(rval_k, rval_kX, pState->res, pState->boxlen, pState->a, pState->om);
        free(rval_kX);
    }
    else if(strcasecmp(name, "L")==0) {
        // angular momentum
        int direction=-1;
        inf >> direction;

        cerr << "Angmom centre is " <<pState->x0 << " " <<pState->y0 << " " << pState->z0 << endl;

        for(long i=0;i<pState->n_part_arr;i++) {
            cen_deriv4_alpha(pState->part_arr[i], direction, rval, pState);
        }
        complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(pState->nPartTotal,sizeof(complex<MyFloat>));
        fft_r(rval_kX, rval, pState->res, 1);
        // The constraint as derived is on the potential. By considering
        // unitarity of FT, we can FT the constraint to get the constraint
        // on the density.
        poiss(rval_k, rval_kX, pState->res, pState->boxlen, pState->a, pState->om);
        free(rval_kX);
    } else {
        cout << "  -> UNKNOWN constraint vector type, returning zeros [this is bad]" << endl;

    }

    if(pState->pField_x!=NULL) {
        cout << " dot in real space = " << std::real(dot(pState->pField_x, rval, pState->nPartTotal)) << endl;
    }

    free(rval);
    return rval_k;

}

void ensureRealDelta(ic_state *pIcs) {
    if(pIcs->pField_x==NULL) {
        pIcs->pField_x = (complex<MyFloat>*)calloc(pIcs->nPartTotal,sizeof(complex<MyFloat>));
        fft_r(pIcs->pField_x,pIcs->pField_k,pIcs->res,-1);
    }
}


#if 0
int main(int argc, char *argv[]){

  if(argc!=2)
  {
      cerr<<"Usage: ./ICgauss paramfile | Output: Pos,Vel,IDs as hdf5 and/or gadget format, power spectrum and used parameters in textfile."<<endl;
      return -1;
  }

    MyFloat Om0, Ol0, zin, sigma8, Boxlength;
    int out, n, gadgetformat;


	MyFloat in_d=1.0;
	complex<MyFloat> d (in_d,0.); //purely real to avoid headaches
    //read parameterfile
    ifstream inf;
      char flag[11][200];
      char value[15][200];
      int seed;
      string incamb;
      string indir;
     string infs=argv[1];
   ofstream parout;

    inf.open(infs.c_str());
      if(inf.is_open()){
	cout<< "Reading parameter file..." << infs.c_str() << endl;
	if(inf.peek()=='%'){inf.ignore(200, '\n');} //ignores first line, which is a comment starting with %

	inf>> flag[0] >>value[0];
	Om0=atof(value[0]);

	    inf>> flag[1] >>value[1];
	    Ol0=atof(value[1]);

	    inf>> flag[2] >>value[2];
	     sigma8=atof(value[2]);

	    inf >> flag[3]>> value[3];
	    Boxlength=atof(value[3]);

	  inf>> flag[4] >>value[4];
	   zin=atof(value[4]);

	  inf>> flag[5] >>value[5];
	   n=atoi(value[5]);

	  inf>> flag[6] >>value[6];
	   out=atoi(value[6]);

	  inf>> flag[7] >>value[7];
	  seed=atoi(value[7]);
	  inf>> flag[8] >>value[8];

	  inf>> flag[9] >>value[9];
	  string indir=value[9];

	  inf>> flag[10] >>value[10];
          gadgetformat=atoi(value[10]);


      }
      else{cerr<< "Noooooo! You fool! That file is not here!"<< endl; return -1;}


     long nPartTotal=(long) (n*n*n);


     gsl_rng * r;
     const gsl_rng_type * T; //generator type variable

     //gsl_rng_env_setup(); //use this instead of gsl_rng_set(r,seed) below for seed from command line

#ifndef DOUBLEPRECISION
T = gsl_rng_ranlxs2; //this is the name of the generator defined in the environmental variable GSL_RNG_TYPE
#else
T = gsl_rng_ranlxs2; //double precision generator: gsl_rng_ranlxd2 //TODO decide which one (but better consistently)
#endif
       r = gsl_rng_alloc (T); //this allocates memory for the generator with type T
       gsl_rng_set(r,seed);


      string base=make_base( value[9], n, Boxlength, zin);
      cout<< "Writing output to "<< (base+ ".*").c_str() <<endl;

      parout.open((base+ ".params").c_str());
      int j;
      //TODO adjust this to write out everything (automatically, instead of setting jmax by hand)
      for(j=0;j<9;j++){parout<< flag[j]<< "	" <<value[j]<<endl;}
      parout.close();

     MyFloat ain=1./(MyFloat)(zin+1);
      int quoppas=600; //max. lines in camb power spectrum file
      int c=7; //for transfer function

      double *inarr=(double*)calloc(quoppas*c,sizeof(double));
      double *kcamb=(double*)calloc(quoppas,sizeof(double));
      double *Tcamb=(double*)calloc(quoppas,sizeof(double));

      if(out!=0 && out!=1 && out!=2){cerr<< "Wrong output format, choose 0 (HDF5), 1 (Gadget) or 2 (both)"<<endl; return -1;}

#ifndef HAVE_HDF5
  if(out!=1){cerr<< "Not compiled with HDF5. Only output=1 is allowed in this case!"<<endl; return -1;}
#endif

      cout<< "Reading transfer file "<< value[8] <<"..."<<endl;

	GetBuffer(inarr, value[8], quoppas*c);

      MyFloat ap=inarr[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

      int quoppa=0;

      for(j=0;j<quoppas;j++){
	if(inarr[c*j]>0){kcamb[j]=MyFloat(inarr[c*j]); Tcamb[j]=MyFloat(inarr[c*j+1])/MyFloat(ap); quoppa+=1;}
	else {continue;}
      }

	free(inarr);

       complex<MyFloat> *rnd=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
       complex<MyFloat> *ft=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

       cout<< "Drawing random numbers..."<< endl;

       long i;
       MyFloat sigma=sqrt((MyFloat)(nPartTotal));
       for(i=0;i<nPartTotal;i++){rnd[i]=gsl_ran_gaussian_ziggurat(r,1.)*sigma;}// cout<< "rnd "<< rnd[i] << endl;}



       gsl_rng_free (r);


       cout<< "First FFT..." <<endl;

       ft=fft_r(ft, rnd, n, 1);

       std::cout<<"Initial chi^2 (white noise, real space) = " << dot(rnd,rnd,nPartTotal)/((MyFloat) nPartTotal) << std::endl;
       free(rnd);

       ft[0]=complex<MyFloat>(0.,0.); //assure mean==0


       std::cout<<"Initial chi^2 (white noise, fourier space) = " << dot(ft,ft,nPartTotal)/((MyFloat) nPartTotal) << std::endl;

      //scale white-noise delta with initial PS
      complex<MyFloat> *ftsc=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      int ix,iy,iz,idx;
      int iix, iiy, iiz;
      int res=n;
      MyFloat kfft;
      MyFloat grwfac;
      MyFloat ns=0.96; //TODO make this an input parameter (long term: run CAMB with input parameters to ensure consistency?)



      //growth factor normalised to 1 today:
      grwfac=D(ain, Om0, Ol0)/D(1., Om0, Ol0);

      cout<< "Growth factor " << grwfac << endl;
      MyFloat sg8;

      sg8=sig(8., kcamb, Tcamb, ns, Boxlength, n, quoppa);
      std::cout <<"Sigma_8 "<< sg8 << std::endl;

      MyFloat kw = 2.*M_PI/(MyFloat)Boxlength;
      MyFloat amp=(sigma8/sg8)*(sigma8/sg8)*grwfac*grwfac; //norm. for sigma8 and linear growth factor
      MyFloat norm=kw*kw*kw/powf(2.*M_PI,3.); //since kw=2pi/L, this is just 1/V_box

      complex<MyFloat> *P=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));


      std::cout<<"Interpolation: kmin: "<< kw <<" Mpc/h, kmax: "<< (MyFloat)( kw*res/2.*sqrt(3.)) <<" Mpc/h"<<std::endl;

       MyFloat norm_amp=norm*amp;
        brute_interpol_new(n, kcamb, Tcamb,  quoppa, kw, ns, norm_amp, ft, ftsc, P);

      //assert(abs(real(norm_iter)) >1e-12); //norm!=0
      //assert(abs(imag(norm_iter)) <1e-12); //because there shouldn't be an imaginary part since we assume d is purely real
      //TODO think about this (it fails more often than it should?)

      cout<<"Transfer applied!"<<endl;
      cout<< "Power spectrum sample: " << P[0] << " " << P[1] <<" " << P[nPartTotal-1] <<endl;

      std::cout<<"Initial chi^2 (white noise, fourier space) = " << chi2(ftsc,P,nPartTotal) << std::endl;






       Grid grid(n);

       // construct the description of the underlying field and its realization
       UnderlyingField<MyFloat> *pField = new UnderlyingField<MyFloat>(P, ftsc, nPartTotal);
       MultiConstrainedField<MyFloat> constr(pField, nPartTotal);

       MyFloat orig_chi2 = std::real(chi2(ftsc,P,nPartTotal));
       std::cout<<"Initial chi^2 (scaled) = " << orig_chi2 << std::endl;

       ic_state ics;
       ics.boxlen = Boxlength;
       ics.dx = Boxlength/n;
       ics.res = n;
       ics.a = ain;
       ics.om=Om0;
       ics.pGrid = &grid;
       ics.nPartTotal=nPartTotal;
       ics.part_arr=NULL;
       ics.n_part_arr=0;
       ics.pField_k = ftsc;
       ics.pField_x = NULL;

       while(!inf.eof()) {
           char command[100];
           command[0]=0;
           inf >> command;

           if(strcasecmp(command,"IDfile")==0) {
               char IDfile[100];
               inf >> IDfile;
               AllocAndGetBuffer_int(IDfile, &ics);
               cout << "New particle array " << IDfile << " loaded." << endl;
           } else if(strcasecmp(command,"append_IDfile")==0) {
               char IDfile[100];
               inf >> IDfile;
               AllocAndGetBuffer_int(IDfile, &ics, true);
           } else if(strcasecmp(command,"select_sphere")==0) {
               float radius;
               inf >> radius;
               selectSphere(&ics,radius);
           }else if(strcasecmp(command,"centre_max")==0) {
               centreDenmax(&ics);
           } else if(strcasecmp(command,"centre_on")==0) {
               long id;
               inf >> id;
               centreOn(&ics, id);
           }
           else if(strcasecmp(command,"order")==0) {
               reorderBuffer(&ics);
           } else if(strcasecmp(command,"truncate")==0) {
               float x;
               inf >> x;
               if(x<0 ) {
                   cerr << "Truncate command takes a fraction between 0 and 1 or a number>=1" << endl;
                   exit(0);
               }
               if(x<1)
                   ics.n_part_arr = ((int)(ics.n_part_arr*x));
               else
                   ics.n_part_arr = ((int)x);
           }
           else if(strcasecmp(command,"calculate")==0) {
               complex<MyFloat> *vec = calcConstraintVector(inf, &ics);
               cout << "    --> calculated value = " << dot(vec, ftsc, nPartTotal) << endl;
               free(vec);
           } else
	   if(strcasecmp(command,"constrain_direction")==0) {
	     // syntax: constrain_direction [and_renormalize] vec_name dir0 dir1 dir2 [renorm_fac]
	     std::stringstream ss;
	     char name[100];
	     bool normalization=false;
	     complex<MyFloat>* vecs[3];
	     complex<MyFloat> vals[3];
	     inf >> name;
	     if(strcasecmp(name,"and_renormalize")==0) {
	       normalization=true;
	       inf >> name;
	     }
	     for(int dir=0; dir<3; dir++) {
	       ss << name << " " ;
	       ss << dir << " " ;
	       vecs[dir] = calcConstraintVector(ss, &ics);
	       vals[dir] = dot(vecs[dir],ftsc,nPartTotal);
	     }

	     complex<MyFloat> norm = std::sqrt(dot(vals,vals,3));
	     cerr << "   Initial values are " << vals[0] << " " << vals[1] << " " <<vals[2] << " -> norm = " << norm << endl;
	     MyFloat direction[3];
	     inf >> direction[0] >> direction[1] >> direction[2];
	     complex<MyFloat> in_norm = std::sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);

	     MyFloat costheta=0;
	     for(int dir=0; dir<3; dir++)
	       costheta+=direction[dir]*std::real(vals[dir])/std::real(norm*in_norm);

	     cerr << "   Between Re original and Re constrained, cos theta = " <<costheta << endl;

	     if(normalization) {
	       MyFloat renorm;
	       inf >> renorm;
	       norm*=renorm;
	     }

	     for(int dir=0; dir<3; dir++)
	       vals[dir]=direction[dir]*norm/in_norm;


	     cerr << "   Constrain values are " << vals[0] << " " << vals[1] << " " << vals[2] << endl;
	     for(int dir=0; dir<3; dir++)
	       constr.add_constraint(vecs[dir],vals[dir],vals[dir]*in_norm/norm);
	   } else
           if(strcasecmp(command,"constrain")==0) {
               bool relative=false;
               if(pField==NULL) {
                   cerr << "Eek! You're trying to add a constraint but the calculation is already done. Move your done command." << endl;
                   exit(0);
               }

               complex<MyFloat> *vec = calcConstraintVector(inf, &ics);

               inf >> command;
               if (strcasecmp(command,"relative")==0) {
                   relative=true;
               } else if (strcasecmp(command,"absolute")!=0) {
                   cerr << "Constraints must state either relative or absolute" << endl;
                   exit(0);
               }

               MyFloat constraint_real;
               inf >> constraint_real;

               std::complex<MyFloat> constraint = constraint_real;
               std::complex<MyFloat> initv = dot(vec, ftsc, nPartTotal);

               if(relative) constraint*=initv;

               cout << "    --> initial value = " << initv << ", constraining to " << constraint << endl;
	           constr.add_constraint(vec, constraint, initv);
            } else if (strcasecmp(command,"cov")==0) {
                constr.print_covariance();
            }else
           if(strcasecmp(command,"done")==0) {
               if(pField==NULL) {
                   cerr << "ERROR - field information is not present. Are there two 'done' commands perhaps?" << endl;
                   exit(0);
               }
               complex<MyFloat> *y1=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
	           constr.prepare();
               constr.get_realization(y1);

               for(long i=0; i<nPartTotal; i++)
                   ftsc[i]=y1[i];

               free(y1);

               cout << "Expected Delta chi^2=" << constr.get_delta_chi2() << endl;

               // clean everything up
               delete pField;
               pField=NULL;

           }


       }

       if(pField!=NULL) {
           cerr << endl << endl << "WHOOPS - you didn't actually calculate the constraints. You need a 'done' command in the paramfile." << endl << endl;
       }


       inf.close();
       MyFloat final_chi2 = std::real(chi2(ftsc,P,nPartTotal));
       std::cout<<"Final chi^2  = " <<  final_chi2 << std::endl;
       std::cout<<"Delta chi^2  = " <<  final_chi2 - orig_chi2 << std::endl;


      cout << endl;




   // ftsc now contains constrained field


   //output power spectrum of constrained field
   powsp_noJing(n, ftsc, (base+ ".ps").c_str(), Boxlength);


   //calculate potential of constrained field
   complex<MyFloat>* potk=(complex<MyFloat>*)calloc(n*n*n,sizeof(complex<MyFloat>));
   potk=poiss(potk, ftsc, n, Boxlength, ain, Om0); //pot in k-space

   complex<MyFloat>* pot=(complex<MyFloat>*)calloc(n*n*n,sizeof(complex<MyFloat>));
   pot=fft_r(pot,potk,res,-1); //pot in real-space



      free(ft);
      free(kcamb);
      free(Tcamb);
      free(P);


      complex<MyFloat>* psift1k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      complex<MyFloat>* psift2k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      complex<MyFloat>* psift3k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      complex<MyFloat>* psift1=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      complex<MyFloat>* psift2=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      complex<MyFloat>* psift3=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

      //complex<MyFloat>* test_arr=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

      //Zeldovich approx.
      for(ix=0; ix<res;ix++){
	for(iy=0;iy<res;iy++){
	  for(iz=0;iz<res;iz++){
	    idx = (ix*res+iy)*(res)+iz;

	      //test_arr[idx].real(MyFloat(idx));
	      if( ix>res/2 ) iix = ix - res; else iix = ix;
	      if( iy>res/2 ) iiy = iy - res; else iiy = iy;
	      if( iz>res/2 ) iiz = iz - res; else iiz = iz;

	      kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);

	      psift1k[idx].real(-ftsc[idx].imag()/(MyFloat)(kfft*kfft)*iix/kw);
	      psift1k[idx].imag(ftsc[idx].real()/(MyFloat)(kfft*kfft)*iix/kw);
	      psift2k[idx].real(-ftsc[idx].imag()/(MyFloat)(kfft*kfft)*iiy/kw);
	      psift2k[idx].imag(ftsc[idx].real()/(MyFloat)(kfft*kfft)*iiy/kw);
	      psift3k[idx].real(-ftsc[idx].imag()/(MyFloat)(kfft*kfft)*iiz/kw);
	      psift3k[idx].imag(ftsc[idx].real()/(MyFloat)(kfft*kfft)*iiz/kw);
	  }
	}
      }



      complex<MyFloat> *delta_real=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
      delta_real=fft_r(delta_real,ftsc,res,-1);

      MyFloat *potr=(MyFloat*)calloc(n*n*n,sizeof(MyFloat));
      for(i=0;i<n*n*n;i++){potr[i]=(MyFloat)(pot[i].real()); pot[i].imag(0.);}

      //optional: save potential
      //string phases_out=(base+ "_phases.hdf5").c_str();
      //save_phases(potk, potr, delta_real, n, phases_out.c_str());

      MyFloat exp=0.;
      for(i=0;i<n*n*n;i++){exp+=(MyFloat)(pot[i].real()*pot[i].real());}
      exp/=(MyFloat)(n*n*n);

      free(pot);
      free(potk);
      free(potr);

      free(delta_real);

      free(ftsc);

       psift1k[0]=complex<MyFloat>(0.,0.);
       psift2k[0]=complex<MyFloat>(0.,0.);
       psift3k[0]=complex<MyFloat>(0.,0.);

       psift1=fft_r(psift1,psift1k,n,-1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
       psift2=fft_r(psift2,psift2k,n,-1); //same
       psift3=fft_r(psift3,psift3k,n,-1); //same

      free(psift1k);
      free(psift2k);
      free(psift3k);



   if (gadgetformat==2){
       SaveGadget2(CreateGadget2Header(nPartTotal, pmass, ain, zin, Boxlength, Om0, Ol0))
     }

    if (gadgetformat==3){

     }

	string dname="Coordinates";
	string dname2="Velocities";
	string dname3="ParticleIDs";


	if(out==0){
#ifdef HAVE_HDF5
	    SaveHDF( (base+ ".hdf5").c_str(), n, header3, Pos1, Pos2, Pos3, Vel1,  Vel2,  Vel3, dname.c_str(), dname2.c_str(), dname3.c_str());
#endif
	}

	else if(out==1){
	    if (gadgetformat==2) SaveGadget2( (base+ "_gadget2.dat").c_str(), n, header2, Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
	    if (gadgetformat==3) SaveGadget3( (base+ "_gadget3.dat").c_str(), n, header3, Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
	}

	else{
	    if (gadgetformat==2) SaveGadget2( (base+ "_gadget2.dat").c_str(), n, header2, Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
	    if (gadgetformat==3) SaveGadget3( (base+ "_gadget3.dat").c_str(), n, header3, Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
#ifdef HAVE_HDF5
	    SaveHDF( (base+ ".hdf5").c_str(), n, header3, Pos1, Pos2, Pos3, Vel1,  Vel2,  Vel3,  dname.c_str(), dname2.c_str(), dname3.c_str());
#endif
	}

	free(Pos1);
	free(Vel1);
	free(Pos2);
	free(Vel2);
	free(Pos3);
	free(Vel3);

       cout<<"Done!"<<endl;

       return 0;
}
#endif
