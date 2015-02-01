#include <string>
#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

template<typename MyFloat>
struct floatinfo {
    static constexpr char const*  name="unknown";
};

template<>
struct floatinfo<double> {
    static constexpr char const*  name="doub";
};

template<>
struct floatinfo<float> {
    static constexpr char const*  name="float";
};

template<typename MyFloat>
class IC {
private:

    MyFloat Om0, Ol0, zin, ain, sigma8, Boxlength;
    int out, n, gadgetformat, seed;

    int quoppa; // eh?

    double *kcamb, *Tcamb;

    long npartTotal;

    string incamb, indir, base;

    complex<MyFloat> *ftsc;

    Grid grid;

public:
    IC() {}

    void setOmegaM0(MyFloat in) {
        Om0=in;
    }

    void setOmegaLambda0(MyFloat in) {
        Ol0=in;
    }

    void setSigma8(MyFloat in) {
        sigma8 = in;
    }

    void setBoxLen(MyFloat in) {
        Boxlength = in;
    }

    void setZ0(MyFloat in) {
        zin = in;
        ain=1./(zin+1.);
    }

    void setn(int in) {
        n = in;
    }

    void setOutputMode(int in) {
        out = in; // the joy of writing this line is substantial

        if(out!=0 && out!=1 && out!=2)
            throw runtime_error("Wrong output format, choose 0 (HDF5), 1 (Gadget) or 2 (both)");

#ifndef HAVE_HDF5
        if(out!=1)
            throw runtime_error("Not compiled with HDF5. Only output=1 is allowed in this case!");
#endif

    }

    void setSeed(int in) {
        seed = in;
    }

    void setCambDat(std::string in) {
        incamb = in;
    }

    void setOutDir(std::string in) {
        indir = in;
    }

    void setGadgetFormat(int in) {
        gadgetformat = in;
        basicRun();
    }

    string make_base( string basename, int res, MyFloat box, MyFloat zin){
      ostringstream result;
      result << basename<<"IC_iter_" << floatinfo<MyFloat>::name << "_z"<<zin<<"_"<<res<<"_L" << box;
      return result.str();
    }

    void readCamb() {

        int quoppas=600; //max. lines in camb power spectrum file
        int c=7; //for transfer function

        double *inarr=(double*)calloc(quoppas*c,sizeof(double));
        kcamb=(double*)calloc(quoppas,sizeof(double));
        Tcamb=(double*)calloc(quoppas,sizeof(double));

        cerr << "Reading transfer file "<< incamb << "..." << endl;
        GetBuffer(inarr, incamb.c_str(), quoppas*c);
        MyFloat ap=inarr[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

        quoppa=0;

        for(int j=0;j<quoppas;j++)
        {
          if(inarr[c*j]>0){kcamb[j]=MyFloat(inarr[c*j]); Tcamb[j]=MyFloat(inarr[c*j+1])/MyFloat(ap); quoppa+=1;}
          else {continue;}
        }

        free(inarr);

    }

    void drawRandom() {
        gsl_rng * r;
        const gsl_rng_type * T; //generator type variable

        T = gsl_rng_ranlxs2; // shouldn't this be gsl_rng_ranlxd2 for MyFloat = double?
        r = gsl_rng_alloc (T); //this allocates memory for the generator with type T
        gsl_rng_set(r,seed);

        complex<MyFloat> *rnd=(complex<MyFloat>*)calloc(npartTotal,sizeof(complex<MyFloat>));
        ftsc=(complex<MyFloat>*)calloc(npartTotal,sizeof(complex<MyFloat>));

        cerr << "Drawing random numbers..."<< endl;

        long i;
        MyFloat sigma=sqrt((MyFloat)(npartTotal));
        for(i=0;i<npartTotal;i++){rnd[i]=gsl_ran_gaussian_ziggurat(r,1.)*sigma;}// cout<< "rnd "<< rnd[i] << endl;}

        gsl_rng_free (r);

        cout<< "First FFT..." <<endl;
        fft_r(ftsc, rnd, n, 1);

        free(rnd);

        ftsc[0]=complex<MyFloat>(0.,0.); //assure mean==0

    }

    void applyPowerSpec() {
        //scale white-noise delta with initial PS

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

        complex<MyFloat> *P=(complex<MyFloat>*)calloc(npartTotal,sizeof(complex<MyFloat>));


        std::cout<<"Interpolation: kmin: "<< kw <<" Mpc/h, kmax: "<< (MyFloat)( kw*res/2.*sqrt(3.)) <<" Mpc/h"<<std::endl;

        MyFloat norm_amp=norm*amp;

        brute_interpol_new(n, kcamb, Tcamb,  quoppa, kw, ns, norm_amp, ftsc, ftsc, P);

        //assert(abs(real(norm_iter)) >1e-12); //norm!=0
        //assert(abs(imag(norm_iter)) <1e-12); //because there shouldn't be an imaginary part since we assume d is purely real
        //TODO think about this (it fails more often than it should?)

        cout<<"Transfer applied!"<<endl;
        cout<< "Power spectrum sample: " << P[0] << " " << P[1] <<" " << P[npartTotal-1] <<endl;

        std::cout<<"Initial chi^2 (white noise, fourier space) = " << chi2(ftsc,P,npartTotal) << std::endl;
    }

    void basicRun() {
        npartTotal = ((long)n*n)*n;

        base = make_base(indir,n,Boxlength,zin);

        readCamb();
        drawRandom();
        applyPowerSpec();

        grid = Grid(n);
    }


    void write() {

    }

/*
    dispatch.add_class_route("IDfile",IC::loadID);
    dispatch.add_class_route("append_IDfile",IC::appendID);
    dispatch.add_class_route("select_sphere",IC::selectSphere);
    dispatch.add_class_route("centre_max",IC::centreDenmax);
    dispatch.add_class_route("centre_on",IC::centreParticle);
    dispatch.add_class_route("order",IC::reorderBuffer);
    dispatch.add_class_route("truncate",IC::truncateBuffer);
    dispatch.add_class_route("calculate",IC::calculate);
    dispatch.add_class_route("constrain",IC::constrain);
    */
};
