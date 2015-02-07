
#ifndef DOUBLEPRECISION     /* default is single-precision */

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
    #ifdef FFTW_TYPE_PREFIX
    #include <drfftw.h>
    #else
    #include <rfftw.h>
    #endif
#endif



template<typename MyFloat>
std::complex<MyFloat> *fft_r(std::complex<MyFloat> *fto, std::complex<MyFloat> *ftin, const int res, const int dir)
{ //works when fto and ftin and output are the same, but overwrites fto for the backwards transform so we have another temporary array there

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
