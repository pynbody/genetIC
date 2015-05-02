#include <string>
#include <tuple>
#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <cassert>
#include <functional>
#include <algorithm>

#include "numpy.hpp"

#define for_each_level(level) for(int level=0; level<2 && n[level]>0; level++)

using namespace std;


template<typename MyFloat>
class IC {
protected:

    MyFloat Om0, Ol0, Ob0, hubble, zin, a, sigma8, ns;

    // Everything about the main grid
    MyFloat boxlen[2], dx[2];
    int n[2]; // number of grid divisions along one axis
    Grid<MyFloat> *pGrid[2]; // the object that helps us relate points on the grid
    complex<MyFloat> *pField_k[2]; // the field in k-space
    complex<MyFloat> *pField_x[2]; // the field in x-space
    complex<MyFloat> *P[2]; // the power spectrum for each k-space point

    complex<MyFloat> *pField_k_0_high; // high-k modes for level 0

    MyFloat x_off[2], y_off[2], z_off[2]; // x,y,z offsets for subgrids

    int zoomfac; // the zoom factor, i.e. the ratio dx/dx[1]
    MyFloat k_cut;   // the cut-off k for low-scale power


    int out, gadgetformat, seed;
    int nCambLines; // eh?
    double *kcamb, *Tcamb;
    long nPartLevel[2];
    long nPartTotal;

    string incamb, indir, inname, base;

    bool whiteNoiseFourier;

    bool prepared;




    std::vector<long> genericParticleArray;
    std::vector<long> levelParticleArray[2];


    std::vector<long> zoomParticleArray;



    MyFloat x0, y0, z0;

    MultiConstrainedField<MyFloat> *pConstrainer;



public:
    IC() {
        pGrid[0]=pGrid[1]=NULL;
        pField_x[0]=pField_x[1]=NULL;
        pField_k[0]=pField_k[1]=NULL;
        pConstrainer=NULL;
        whiteNoiseFourier=false;
        hubble=0.701;   // old default
        Ob0=-1.0;
        ns = 0.96;      // old default
        n[0]=0;
        n[1]=-1; // no subgrid by default
        x_off[0]=y_off[0]=z_off[0]=0;
        x_off[1]=y_off[1]=z_off[1]=0;
        prepared = false;
    }

    ~IC() {
        if(pGrid[0]!=NULL) delete pGrid[0];
        if(pGrid[1]!=NULL) delete pGrid[1];
    }

    void setOmegaM0(MyFloat in) {
        Om0=in;
    }

    void setOmegaB0(MyFloat in) {
        Ob0=in;
    }

    void setOmegaLambda0(MyFloat in) {
        Ol0=in;
    }

    void setHubble(MyFloat in) {
        hubble=in;
    }

    void setSigma8(MyFloat in) {
        sigma8 = in;
    }

    void setBoxLen(MyFloat in) {
        boxlen[0] = in;
        dx[0] = boxlen[0]/n[0];
    }

    void setZ0(MyFloat in) {
        zin = in;
        a=1./(zin+1.);
    }

    void setn(int in) {
        n[0] = in;
        nPartLevel[0] = ((long)n[0]*n[0])*n[0];
        dx[0] = boxlen[0]/n[0];
    }

    void setns(MyFloat in) {
        ns = in;
    }

    void setn2(int in) {
        n[1] = in;
        nPartLevel[1] = ((long)n[1]*n[1])*n[1];
    }

    void setZoom(int in) {
        // open a subgrid which is the specified factor smaller than
        // the parent grid
        boxlen[1] = boxlen[0]/in;
        zoomfac = in;
    }

    void setZoomParticles(string fname) {
        if(n[1]==0)
            throw(std::runtime_error("Set n2 before specifying the zoom particles"));

        if(boxlen[1]==0)
            throw(std::runtime_error("Set the zoom factor before specifying the zoom particles"));

        AllocAndGetBuffer_int(fname.c_str());
        zoomParticleArray = genericParticleArray;

        genericParticleArray.clear();

        // Not strictly necessary for this function, but we rely on having the
        // particle IDs sorted later on (when writing)
        std::sort(zoomParticleArray.begin(), zoomParticleArray.end());



        // find boundaries
        int x0, x1, y0, y1, z0, z1;
        int x,y,z;

        x0=y0=z0=n[0];
        x1=y1=z1=0;
        Grid<MyFloat> g(n[0]);

        // TO DO: wrap the box sensibly

        for(long i=0; i<zoomParticleArray.size(); i++) {
            g.get_coordinates(zoomParticleArray[i],x,y,z);
            if(x<x0) x0=x;
            if(y<y0) y0=y;
            if(z<z0) z0=z;
            if(x>x1) x1=x;
            if(y>y1) y1=y;
            if(z>z1) z1=z;
            // cerr << i << " " << zoomParticleArray[i] << " " << x << " " << y << " " << z << endl;
        }

        // Now see if the zoom the user chose is OK
        int n_user = n[0]/zoomfac;
        if((x1-x0)>n_user || (y1-y0)>n_user || (z1-z0)>n_user) {
            cerr << x1 << " " << x0 << " " << n_user << endl;
            throw(std::runtime_error("Zoom particles do not fit in specified sub-box. Decrease zoom, or choose different particles. (NB wrapping not yet implemented)"));
        }
        // At this point we know things fit. All we need to do is choose
        // the correct offset to get the particles near the centre of the
        // zoom box.

        // Here is the bottom left of the box:
        x=(x0+x1)/2-n[0]/(2*zoomfac);
        y=(y0+y1)/2-n[0]/(2*zoomfac);
        z=(z0+z1)/2-n[0]/(2*zoomfac);

        if(x<0) x=0;
        if(y<0) y=0;
        if(z<0) z=0;
        if(x>n[0]-n[0]/zoomfac) x = n[0]-n[0]/zoomfac;
        if(y>n[0]-n[0]/zoomfac) y = n[0]-n[0]/zoomfac;
        if(z>n[0]-n[0]/zoomfac) z = n[0]-n[0]/zoomfac;


        x_off[1] = x_off[0] + x*dx[0];
        y_off[1] = y_off[0] + y*dx[0];
        z_off[1] = z_off[0] + z*dx[0];



        initZoom();

    }

    void initZoom() {
        zoomfac = (boxlen[0]/n[0])/(boxlen[1]/n[1]);
        dx[1] = boxlen[1]/n[1];
        nPartTotal = nPartLevel[0]+(zoomfac*zoomfac*zoomfac-1)*zoomParticleArray.size();
        cout << "Initialized a zoom region:" << endl;
        cout << "  Subbox length   = " << boxlen[1] << " Mpc/h" << endl;
        cout << "  n[1]            = " << n[1] << endl;
        cout << "  dx[1]           = " << dx[1] << endl;
        cout << "  Zoom factor     = " << zoomfac << endl;
        cout << "  Low-left corner = " << x_off[1] << ", " << y_off[1] << ", " << z_off[1] << endl;
        cout << "  Total particles = " << nPartTotal << endl;

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

    void setSeedFourier(int in) {
        seed = in;
        whiteNoiseFourier = true;
    }

    void setCambDat(std::string in) {
        incamb = in;
    }

    void setOutDir(std::string in) {
        indir = in;
    }

    void setOutName(std::string in) {
        inname=in;
    }

    void setGadgetFormat(int in) {
        gadgetformat = in;
        if(!prepared)
            prepare(); //compatibility with old paramfiles
    }

    string make_base( string basename, int level=0){
        ostringstream nult;
        if(inname.size()==0) {
            nult << basename<<"IC_iter_" << floatinfo<MyFloat>::name << "_z"<<zin<<"_"<<n[level]<<"_L" << boxlen[level];

        } else {
            if (level==0)
                nult << basename << "/" << inname;
            else
                nult << basename << "/" << inname << "_" << level;
        }
        return nult.str();
    }

    void readCamb() {

        const int maxlines=600; //max. lines in camb power spectrum file
        const int c=7; //for transfer function
        int j;
        double *inarr=(double*)calloc(maxlines*c,sizeof(double));
        kcamb=(double*)calloc(maxlines,sizeof(double));
        Tcamb=(double*)calloc(maxlines,sizeof(double));

        cerr << "Reading transfer file "<< incamb << "..." << endl;
        GetBuffer(inarr, incamb.c_str(), maxlines*c);
        MyFloat ap=inarr[1]; //to normalise CAMB transfer function so T(0)= 1, doesn't matter if we normalise here in terms of accuracy, but feels more natural

        nCambLines=0;

        for(j=0;j<maxlines;j++)
        {
          if(inarr[c*j]>0){kcamb[j]=MyFloat(inarr[c*j]); Tcamb[j]=MyFloat(inarr[c*j+1])/MyFloat(ap); nCambLines+=1;}
          else {continue;}
        }

        // extend high-k range using power law

        MyFloat gradient = log(Tcamb[nCambLines-1]/Tcamb[nCambLines-2])/log(kcamb[nCambLines-1]/kcamb[nCambLines-2]);
        cerr << "Extending CAMB transfer using powerlaw " << gradient << " from " << Tcamb[nCambLines-1] << endl;

        for(j=nCambLines;kcamb[j-1]<100;j++)
        {
            kcamb[j] = kcamb[j-1]+1.0;
            Tcamb[j] = exp(log(Tcamb[nCambLines-1]) + gradient * (log(kcamb[j]/kcamb[nCambLines-1])));
            // cerr << kcamb[j-1] << " " << Tcamb[j-1] << endl;
        }
        nCambLines=j-1;

        free(inarr);

    }


    static long kgridToIndex(int k1, int k2, int k3, int n)  {
        long ii,jj,ll,idk;
        if(k1<0) ii=k1+n; else ii=k1;
        if(k2<0) jj=k2+n; else jj=k2;
        if(k3<0) ll=k3+n; else ll=k3;
        idk=(ii*(n)+jj)*(n)+ll;
        return idk;
    }

    static void drawOneFourierMode(gsl_rng *r, int k1, int k2, int k3,
                                  MyFloat norm, int n, complex<MyFloat> *pField_k) {
        long id_k, id_negk;
        id_k = kgridToIndex(k1,k2,k3,n);
        id_negk = kgridToIndex(-k1,-k2,-k3,n);
        pField_k[id_k]=std::complex<MyFloat>(norm*gsl_ran_gaussian_ziggurat(r,1.),norm*gsl_ran_gaussian_ziggurat(r,1.));

        // reality condition:
        pField_k[id_negk]=std::conj(pField_k[id_k]);
    }


    static void drawRandomForSpecifiedGrid(int n, complex<MyFloat> *& pField_x,
                                           complex<MyFloat> *& pField_k,
                                           gsl_rng * r)
    {
        long nPartTotal = n;
        nPartTotal*=n*n;

        pField_k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        MyFloat sigma=sqrt((MyFloat)(nPartTotal));

        complex<MyFloat> *rnd=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

        cerr << "Drawing random numbers..."<< endl;

        long i;

        for(i=0;i<nPartTotal;i++){rnd[i]=gsl_ran_gaussian_ziggurat(r,1.)*sigma;}// cout<< "rnd "<< rnd[i] << endl;}



        cout<< "First FFT..." <<endl;
        fft_r(pField_k, rnd, n, 1);

        free(rnd);

        pField_k[0]=complex<MyFloat>(0.,0.); //assure mean==0



    }

    static void drawRandomForSpecifiedGridFourier(int n,
                                           complex<MyFloat> *& pField_k,
                                           gsl_rng * r) {

       long nPartTotal = n;
       nPartTotal*=n*n;

       pField_k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
       MyFloat sigma=sqrt((MyFloat)(nPartTotal));


       // do it in fourier space, in order of increasing |k|, so that
       // resolution can be scaled by factors of 2 and we still get
       // the 'same' field

       cerr << "Drawing random numbers in fourier space..."<< endl;
       MyFloat kk=0.;
       int ks,k1,k2,k3;

       sigma/=sqrt(2.0);

       // Do it in square k-shells
       for(ks=0; ks<n/2;ks++) {
           for(k1=-ks; k1<ks; k1++) {
               for(k2=-ks; k2<ks; k2++) {
                   drawOneFourierMode(r,ks,k1,k2,sigma,n,pField_k);
                   drawOneFourierMode(r,k1,ks,k2,sigma,n,pField_k);
                   drawOneFourierMode(r,k1,k2,ks,sigma,n,pField_k);
               }
           }
        }

    }

    void drawRandom() {
        gsl_rng * r;
        const gsl_rng_type * T; //generator type variable

        T = gsl_rng_ranlxs2; // shouldn't this be gsl_rng_ranlxd2 for MyFloat = double?
        r = gsl_rng_alloc (T); //this allocates memory for the generator with type T
        gsl_rng_set(r,seed);

        for_each_level(level) {
            cerr << "Generate random noise for level " << level <<endl;
            if(whiteNoiseFourier)
                drawRandomForSpecifiedGridFourier(n[level], pField_k[level],r);
            else
                drawRandomForSpecifiedGrid(n[level], pField_x[level], pField_k[level],r);
        }

        gsl_rng_free (r);
    }

    void zeroLevel(int level) {
        cerr << "*** Warning: your script calls zeroLevel("<<level <<"). This is intended for testing purposes only!" << endl;
        if(level==-1) {
            if(pField_k_0_high) memset(pField_k_0_high,0,sizeof(MyFloat)*2*nPartLevel[0]);
        } else {
            if(pField_x[level]!=NULL) {
                memset(pField_x[level],0,sizeof(MyFloat)*2*nPartLevel[level]);
            }
            if(pField_k[level]!=NULL) {
                memset(pField_k[level],0,sizeof(MyFloat)*2*nPartLevel[level]);
            }
        }
    }

    template<typename T>
    void interpolateLevel(int level, T* pField_x_l, T* pField_x_p) {

        if(level<=0)
            throw std::runtime_error("Trying to interpolate onto the top-level grid");

        int n_l = n[level];
        MyFloat dx_l = dx[level];
        MyFloat dx_p = dx[level-1];

        Grid<MyFloat> *pGrid_l = pGrid[level];
        Grid<MyFloat> *pGrid_p = pGrid[level-1];
        MyFloat x,y,z; // coordinates
        int x_p_0, y_p_0, z_p_0, x_p_1, y_p_1, z_p_1; // parent grid locations
        MyFloat xw0,yw0,zw0; // weights for lower-left parent grid contribution
        MyFloat xw1,yw1,zw1; // weights for upper-right parent grid contribution

        // get offsets of bottom-left of grids
        MyFloat x_off_l = x_off[level], y_off_l = y_off[level], z_off_l = z_off[level];
        MyFloat x_off_p = x_off[level-1], y_off_p = y_off[level-1], z_off_p = z_off[level-1];



        for(int x_l=0; x_l<n_l; x_l++) {
            for(int y_l=0; y_l<n_l; y_l++){
                for(int z_l=0; z_l<n_l; z_l++) {
                    // find centre point of current cell:
                    x = ((MyFloat) x_l+0.5)*dx_l+x_off_l;
                    y = ((MyFloat) y_l+0.5)*dx_l+y_off_l;
                    z = ((MyFloat) z_l+0.5)*dx_l+z_off_l;

                    // grid coordinates of parent cell starting to bottom-left
                    // of our current point
                    x_p_0 = (int) floor(((x-x_off_p)/dx_p - 0.5));
                    y_p_0 = (int) floor(((y-y_off_p)/dx_p - 0.5));
                    z_p_0 = (int) floor(((z-z_off_p)/dx_p - 0.5));

                    // grid coordinates of top-right
                    x_p_1 = x_p_0+1;
                    y_p_1 = y_p_0+1;
                    z_p_1 = z_p_0+1;


                    // weights, which are the distance to the centre point of the
                    // upper-right cell, in grid units (-> maximum 1)

                    xw0 = ((MyFloat) x_p_1 + 0.5) - ((x-x_off_p)/dx_p);
                    yw0 = ((MyFloat) y_p_1 + 0.5) - ((y-y_off_p)/dx_p);
                    zw0 = ((MyFloat) z_p_1 + 0.5) - ((z-z_off_p)/dx_p);

                    xw1 = 1.-xw0;
                    yw1 = 1.-yw0;
                    zw1 = 1.-zw0;

                    assert(xw0<=1.0 && xw0>=0.0);

                    /*
                    cerr << "  " << x_l << " " << y_l << " " << z_l << endl;
                    cerr << "->" << x_p_0 << " " << y_p_0 << " " << z_p_0 << endl;
                    cerr << " w" << xw0 << " " << yw0 << " " << zw0;
                    */

                    long ind_l = pGrid_l->get_index(x_l,y_l,z_l);

                    pField_x_l[ind_l]+= xw0*yw0*zw1*pField_x_p[pGrid_p->get_index(x_p_0,y_p_0,z_p_1)] +
                                        xw1*yw0*zw1*pField_x_p[pGrid_p->get_index(x_p_1,y_p_0,z_p_1)] +
                                        xw0*yw1*zw1*pField_x_p[pGrid_p->get_index(x_p_0,y_p_1,z_p_1)] +
                                        xw1*yw1*zw1*pField_x_p[pGrid_p->get_index(x_p_1,y_p_1,z_p_1)] +
                                        xw0*yw0*zw0*pField_x_p[pGrid_p->get_index(x_p_0,y_p_0,z_p_0)] +
                                        xw1*yw0*zw0*pField_x_p[pGrid_p->get_index(x_p_1,y_p_0,z_p_0)] +
                                        xw0*yw1*zw0*pField_x_p[pGrid_p->get_index(x_p_0,y_p_1,z_p_0)] +
                                        xw1*yw1*zw0*pField_x_p[pGrid_p->get_index(x_p_1,y_p_1,z_p_0)];


                }
            }
        }

    }

    void interpolateLevel(int level) {
        operateInRealSpace(level);
        operateInRealSpace(level-1);
        complex<MyFloat> *pField_x_l = pField_x[level];
        complex<MyFloat> *pField_x_p = pField_x[level-1];
        interpolateLevel(level,pField_x_l,pField_x_p);
    }


    /*
    Much of the griding code is generic and could support 'nested' zoom
    regions. However the following parts are specific to a 2-level system
    */
    MyFloat filter_low(MyFloat k) {
        // Return filter for low-res part of field
        MyFloat k_cut = ((MyFloat) n[0]) * 0.3 * 2.*M_PI/boxlen[0];
        MyFloat T = k_cut/10; // arbitrary
        return 1./(1.+exp((k-k_cut)/T));
    }

    MyFloat filter_high(MyFloat k) {
        // Return filter for high-res part of field
        return 1.-filter_low(k);
    }

    std::function<MyFloat(MyFloat)> C_filter_low() {
        // Return filter for low-res part of covariance
        return [this](MyFloat k){
            MyFloat p = this->filter_low(k);
            return p*p;
        };
    }

    std::function<MyFloat(MyFloat)> C_filter_lowhigh() {
        // Return filter for high-k contribution to low-res part of covariance
        // NB this is sampled from the same original low-res field so the sqrt of
        // the filter needs to sum to one, not the filter itself
        return [this](MyFloat k){
            MyFloat p= this->filter_high(k);
            return p*p;
        };
    }

    std::function<MyFloat(MyFloat)> C_filter_high() {
        // Return filter for high-res part of covariance
        return [this](MyFloat k){
            MyFloat p = this->filter_low(k);
            return MyFloat(1.0)-p*p;
        };
    }

    std::function<MyFloat(MyFloat)> C_filter_default() {
        return
            [](MyFloat k) {
                return MyFloat(1.0);
            };
    }


    void splitLevel0() {
        operateInFourierSpace(0);
        pField_k_0_high = (complex<MyFloat>*)calloc(this->nPartLevel[0],sizeof(complex<MyFloat>));
        memcpy(pField_k_0_high,pField_k[0],sizeof(complex<MyFloat>)*this->nPartLevel[0]);
        cerr << "Apply transfer function for high-k part of level 0 field" << endl;
        applyPowerSpecForLevel(0,true);
    }

    void recombineLevel0() {

        operateInFourierSpace(0);
        for(long i=0; i<this->nPartLevel[0]; i++)
            pField_k[0][i]+=pField_k_0_high[i];
        free(pField_k_0_high);

    }

    void applyPowerSpecForLevel(int level, bool high_k=false) {
        //scale white-noise delta with initial PS

        operateInFourierSpace(level);

        MyFloat grwfac;

        //growth factor normalised to 1 today:
        grwfac=D(a, Om0, Ol0)/D(1., Om0, Ol0);

        cout<< "Growth factor " << grwfac << endl;
        MyFloat sg8;

        sg8=sig(8., kcamb, Tcamb, ns, boxlen[level], n[level], nCambLines);
        std::cout <<"Sigma_8 "<< sg8 << std::endl;

        MyFloat kw = 2.*M_PI/(MyFloat)boxlen[level];

        // TODO: link sigma8 from the top level down

        MyFloat amp=(sigma8/sg8)*(sigma8/sg8)*grwfac*grwfac; //norm. for sigma8 and linear growth factor
        MyFloat norm=kw*kw*kw/powf(2.*M_PI,3.); //since kw=2pi/L, this is just 1/V_box

        P[level]=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));

        std::function<MyFloat(MyFloat)> filter;

        if(level==0 && n[1]>0)
        {
            if(high_k) {
                cerr << "Select C_filter_lowhigh" << endl;
                filter = C_filter_lowhigh();
            } else {
                cerr << "Select C_filter_low" << endl;
                filter = C_filter_low();
            }
        }

        else if(level==1)
        {
            cerr << "Select C_filter_high" << endl;
            filter = C_filter_high();
        }
        else
        {
            filter = C_filter_default();
            cerr << "No filter" << endl;
        }

        std::cout<<"Interpolation: kmin: "<< kw <<" Mpc/h, kmax: "<< (MyFloat)( kw*n[level]/2.*sqrt(3.)) <<" Mpc/h"<<std::endl;

        MyFloat norm_amp=norm*amp;

        complex<MyFloat> *pField_k_this = pField_k[level];
        if(high_k)
            pField_k_this = pField_k_0_high;

        brute_interpol_new(n[level], kcamb, Tcamb,  nCambLines, kw, ns, norm_amp, pField_k_this, pField_k_this, P[level], filter);

        //assert(abs(real(norm_iter)) >1e-12); //norm!=0
        //assert(abs(imag(norm_iter)) <1e-12); //because there shouldn't be an imaginary part since we assume d is purely real
        //TODO think about this (it fails more often than it should?)

        cout<<"Transfer applied!"<<endl;
        // cout<< "Power spectrum sample: " << P[level][0] << " " << P[level][1] <<" " << P[level][nPartTotal-1] <<endl;

        // std::cout<<"Initial chi^2 (white noise, fourier space) = " << chi2(pField_k[level],P[level],nPartTotal) << std::endl;



    }

    void applyPowerSpec() {
        for_each_level(level) {
            cerr << "Apply transfer function for level " << level <<endl;
            applyPowerSpecForLevel(level);
        }
    }

    void dumpGrid(int level=0) {
        ensureRealDelta(level);

        int shape[3] = {n[level], n[level], n[level]};
        int fortran_order = 0;

        ostringstream filename;
        filename << indir << "/grid-" << level << ".npy";

        const int dim[3] = { n[level],n[level],n[level] };
        numpy::SaveArrayAsNumpy(filename.str(), true, 3, dim, pField_x[level] );


        filename.str("");
        filename << indir << "/grid-info-" << level << ".txt";

        ofstream ifile;
        ifile.open(filename.str());
        cerr << "Writing to " << filename.str() << endl;

        ifile << x_off[level] << " " << y_off[level] << " " << z_off[level] << " " << boxlen[level] << endl;
        ifile << "The line above contains information about grid level " << level << endl;
        ifile << "It gives the x-offset, y-offset and z-offset of the low-left corner and also the box length" << endl;
        ifile.close();
    }

    void dumpPS(int level=0) {
        operateInFourierSpace(level);
        powsp_noJing(n[level], pField_k[level], (base+"_"+((char)(level+'0'))+".ps").c_str(), boxlen[level]);
    }

    std::tuple<MyFloat*, MyFloat*, MyFloat*, MyFloat*, MyFloat*, MyFloat*> zeldovich(int level=0, long nPartAllocate=-1) {
        //Zeldovich approx.

        if(nPartAllocate==-1)
            nPartAllocate = nPartLevel[level];

        operateInFourierSpace(level);

        complex<MyFloat>* psift1k=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));
        complex<MyFloat>* psift2k=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));
        complex<MyFloat>* psift3k=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));

        int iix, iiy, iiz;
        long idx;
        MyFloat kfft;
        MyFloat kw = 2.*M_PI/(MyFloat)boxlen[level];

        for(int ix=0; ix<n[level];ix++){
            for(int iy=0;iy<n[level];iy++){
                for(int iz=0;iz<n[level];iz++){


                    idx = (ix*n[level]+iy)*(n[level])+iz;

                    if( ix>n[level]/2 ) iix = ix - n[level]; else iix = ix;
                    if( iy>n[level]/2 ) iiy = iy - n[level]; else iiy = iy;
                    if( iz>n[level]/2 ) iiz = iz - n[level]; else iiz = iz;

                    kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);

                    psift1k[idx].real(-pField_k[level][idx].imag()/(MyFloat)(kfft*kfft)*iix/kw);
                    psift1k[idx].imag(pField_k[level][idx].real()/(MyFloat)(kfft*kfft)*iix/kw);
                    psift2k[idx].real(-pField_k[level][idx].imag()/(MyFloat)(kfft*kfft)*iiy/kw);
                    psift2k[idx].imag(pField_k[level][idx].real()/(MyFloat)(kfft*kfft)*iiy/kw);
                    psift3k[idx].real(-pField_k[level][idx].imag()/(MyFloat)(kfft*kfft)*iiz/kw);
                    psift3k[idx].imag(pField_k[level][idx].real()/(MyFloat)(kfft*kfft)*iiz/kw);
                }
            }
        }


        psift1k[0]=complex<MyFloat>(0.,0.);
        psift2k[0]=complex<MyFloat>(0.,0.);
        psift3k[0]=complex<MyFloat>(0.,0.);

        complex<MyFloat>* psift1=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));
        complex<MyFloat>* psift2=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));
        complex<MyFloat>* psift3=(complex<MyFloat>*)calloc(nPartLevel[level],sizeof(complex<MyFloat>));

        psift1=fft_r(psift1,psift1k,n[level],-1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
        psift2=fft_r(psift2,psift2k,n[level],-1); //same
        psift3=fft_r(psift3,psift3k,n[level],-1); //same

        free(psift1k);
        free(psift2k);
        free(psift3k);

        MyFloat gr=boxlen[level]/(MyFloat)n[level];
        cout<< "Grid cell size: "<< gr <<" Mpc/h"<<endl;

        MyFloat *Vel1=(MyFloat*)calloc(nPartAllocate,sizeof(MyFloat));
        MyFloat *Vel2=(MyFloat*)calloc(nPartAllocate,sizeof(MyFloat));
        MyFloat *Vel3=(MyFloat*)calloc(nPartAllocate,sizeof(MyFloat));
        MyFloat *Pos1=(MyFloat*)calloc(nPartAllocate,sizeof(MyFloat));
        MyFloat *Pos2=(MyFloat*)calloc(nPartAllocate,sizeof(MyFloat));
        MyFloat *Pos3=(MyFloat*)calloc(nPartAllocate,sizeof(MyFloat));

        MyFloat hfac=1.*100.*sqrt(Om0/a/a/a+Ol0)*sqrt(a);
        //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
        //TODO: hardcoded value of f=1 is inaccurate, but fomega currently gives wrong nults


        MyFloat mean1=0.,mean2=0.,mean3=0.;
        MyFloat Mean1=0., Mean2=0., Mean3=0.;
        cout<< "Applying ZA & PBC... "<<endl;
        //apply ZA:
        for(int ix=0;ix<n[level];ix++) {
            for(int iy=0;iy<n[level];iy++) {
                for(int iz=0;iz<n[level];iz++) {

                    idx = (ix*n[level]+iy)*(n[level])+iz;

                    Vel1[idx] = psift1[idx].real()*hfac; //physical units
                    Vel2[idx] = psift2[idx].real()*hfac;
                    Vel3[idx] = psift3[idx].real()*hfac;

                    // position OFFSET in physical coordinates
                    // NB must work in offsets to support interpolation below
                    Pos1[idx] = psift1[idx].real();
                    Pos2[idx] = psift2[idx].real();
                    Pos3[idx] = psift3[idx].real();


                    /*
                    // grid now added elsewhere...

                    Pos1[idx] = fmod(Pos1[idx],boxlen[0]);
                    if(Pos1[idx]<0) Pos1[idx]+=boxlen[0];
                    Pos2[idx] = fmod(Pos2[idx],boxlen[0]);
                    if(Pos2[idx]<0) Pos2[idx]+=boxlen[0];
                    Pos3[idx] = fmod(Pos3[idx],boxlen[0]);
                    if(Pos3[idx]<0) Pos3[idx]+=boxlen[0];


                    Mean1+=Pos1[idx];
                    Mean2+=Pos2[idx];
                    Mean3+=Pos3[idx];
                    */


                }
            }
        }


        free(psift1);
        free(psift2);
        free(psift3);

        return make_tuple(Pos1,Pos2,Pos3,Vel1,Vel2,Vel3);

    }



    void writeLevel(int level=0) {

        MyFloat *Pos1, *Pos2, *Pos3, *Vel1, *Vel2, *Vel3;

        tie(Pos1,Pos2,Pos3,Vel1,Vel2,Vel3) = zeldovich(level);

        pGrid[level]->add_grid(Pos1,Pos2,Pos3,boxlen[0]); // last arg forces wrapping at BASE level always

        MyFloat pmass=27.78*Om0*powf(boxlen[level]/(MyFloat)(n[level]),3.0); // in 10^10 M_sol
        cout<< "Particle mass: " <<pmass <<" [10^10 M_sun]"<<endl;

        if (gadgetformat==2){
            SaveGadget2( (base+ "_gadget2.dat").c_str(), nPartLevel[level],
            CreateGadget2Header<MyFloat>(nPartLevel[level], pmass, a, zin, boxlen[level], Om0, Ol0, hubble),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
        }
        else if (gadgetformat==3){
            SaveGadget3( (base+ "_gadget3.dat").c_str(), nPartLevel[level],
            CreateGadget3Header<MyFloat>(nPartLevel[level], pmass, a, zin, boxlen[level], Om0, Ol0, hubble),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
        }
        else if (gadgetformat==4) {
            SaveTipsy( (base+".tipsy").c_str(), nPartLevel[level], Pos1, Vel1, Pos2, Vel2, Pos3, Vel3,
		       boxlen[level],  Om0,  Ol0,  hubble,  a, pmass, (MyFloat*)NULL, -1.0, 1.0/(140*n[level]));
        }

        free(Pos1);
        free(Vel1);
        free(Pos2);
        free(Vel2);
        free(Pos3);
        free(Vel3);
    }

    ///////////////////////////////////////////////
    // ZOOM particle management
    ///////////////////////////////////////////////

    void deleteParticles(std::vector<MyFloat*> A) {

        // delete the particles in the 'zoom' region.
        // This involves just moving everything backwards
        // in steps.

        long i_zoom=0;
        long write_ptcl=0;
        long next_zoom=zoomParticleArray[0];
        long next_zoom_i=0;

        for(long read_ptcl=0; read_ptcl<nPartLevel[0]; read_ptcl++) {
            if(read_ptcl!=next_zoom) {
                if(read_ptcl!=write_ptcl) {
                    for(auto ar=A.begin(); ar!=A.end(); ++ar) {
                        (*ar)[write_ptcl]=(*ar)[read_ptcl];
                    }
                }
                write_ptcl++;
            } else {
                // we've just hit a 'zoom' particle. Skip over it.
                // Keep track of the next zoom particle
                if(next_zoom_i+1<zoomParticleArray.size()) {
                    next_zoom_i++;
                    next_zoom = zoomParticleArray[next_zoom_i];
                } else {
                    // no more zoom particles
                    next_zoom =-1;
                }
            }
        }
    }

    std::vector<long> mapid(long id0, int level0=0, int level1=1) {
        // Finds all the particles at level1 corresponding to the single
        // particle at level0
        std::tie(x0,y0,z0) = pGrid[level0]->get_centroid_location(id0);
        return pGrid[level1]->get_ids_in_cube(x0,y0,z0,dx[level0]);
    }

    void insertParticles(std::vector<MyFloat*> A, std::vector<MyFloat*> B) {
        // insert the particles from the zoom region


        // the last 'low-res' particle is just before this address:
        long i_write = nPartLevel[0]-zoomParticleArray.size();

        int zoomfac3=zoomfac*zoomfac*zoomfac; // used for testing only

        for(long i=0; i<zoomParticleArray.size(); i++) {
            // get the list of zoomed particles corresponding to this one
            std::vector<long> hr_particles = mapid(zoomParticleArray[i]);
            assert(hr_particles.size()==zoomfac3);
            for(auto i=hr_particles.begin(); i!=hr_particles.end(); i++) {
                for(auto ar_to=A.begin(), ar_from=B.begin();
                    ar_to!=A.end() && ar_from!=B.end();)
                {
                    (*ar_to)[i_write] = (*ar_from)[(*i)];
                    ++ar_from;
                    ++ar_to;
                }
                ++i_write;
            }
        }
        assert(i_write==nPartTotal);
    }

    void interpretParticleList() {

        levelParticleArray[0].clear();
        levelParticleArray[1].clear();

        if(zoomParticleArray.size()==0) {
            levelParticleArray[0] = genericParticleArray;
        } else {

            long i0 = nPartLevel[0]-zoomParticleArray.size();
            long lr_particle;
            int n_hr_per_lr = zoomfac*zoomfac*zoomfac;

            for(auto i=genericParticleArray.begin(); i!=genericParticleArray.end(); ++i) {
                if((*i)<i0)
                    throw std::runtime_error("Constraining particle is in low-res region - not permitted");

                if(((*i)-i0)/n_hr_per_lr>=zoomParticleArray.size())
                    throw std::runtime_error("Particle ID out of range");

                // find the low-res particle. Note that we push it onto the list
                // even if it's already on there to get the weighting right and
                // prevent the need for a costly search
                lr_particle = zoomParticleArray[((*i)-i0)/n_hr_per_lr];

                levelParticleArray[0].push_back(lr_particle);

                // get all the HR particles
                std::vector<long> hr_particles = mapid(lr_particle);
                assert(hr_particles.size()==n_hr_per_lr);

                // work out which of these this particle must be and push it onto
                // the HR list
                int offset = ((*i)-i0)%n_hr_per_lr;
                levelParticleArray[1].push_back(hr_particles[offset]);

            }

        }
    }



    void write() {
        if(n[1]<=0) {
            cerr << "***** WRITE - no zoom *****" << endl;
            writeLevel(0);
        } else {
            cerr << "***** WRITE - zoom, total particles = "<< nPartTotal << " *****" << endl;

            MyFloat *Pos1, *Pos2, *Pos3, *Vel1, *Vel2, *Vel3, *Mass;
            MyFloat *Pos1z, *Pos2z, *Pos3z, *Vel1z, *Vel2z, *Vel3z;

            tie(Pos1,Pos2,Pos3,Vel1,Vel2,Vel3) = zeldovich(0);

            MyFloat pmass1=27.78*Om0*powf(boxlen[0]/(MyFloat)(n[0]),3.0); // in 10^10 M_sol
            MyFloat pmass2 = 27.78*Om0*powf(boxlen[1]/(MyFloat)(n[1]),3.0); // in 10^10 M_sol

            cerr<< "Main particle mass: " <<pmass1 <<" [10^10 M_sun]"<<endl;
            cerr<< "Zoom particle mass: " <<pmass2 <<" [10^10 M_sun]"<<endl;

            tie(Pos1z,Pos2z,Pos3z,Vel1z,Vel2z,Vel3z) = zeldovich(1);

            // Interpolate the low-frequency information from level 0:
            interpolateLevel(1,Pos1z,Pos1);
            interpolateLevel(1,Pos2z,Pos2);
            interpolateLevel(1,Pos3z,Pos3);
            interpolateLevel(1,Vel1z,Vel1);
            interpolateLevel(1,Vel2z,Vel2);
            interpolateLevel(1,Vel3z,Vel3);

            // This isn't strictly necessary but it can make for useful output:
            // interpolate the density field
            interpolateLevel(1);

            // Now re-do level 0, but include the high-k modes which have
            // earlier been filtered out
            free(Pos1); free(Pos2); free(Pos3); free(Vel1); free(Vel2); free(Vel3);

            recombineLevel0();

            tie(Pos1,Pos2,Pos3,Vel1,Vel2,Vel3) = zeldovich(0, nPartTotal);


            // add the grid offsets to each position:
            pGrid[0]->add_grid(Pos1,Pos2,Pos3,boxlen[0]);
            pGrid[1]->add_grid(Pos1z,Pos2z,Pos3z,boxlen[0]);

            // now we go through and populate the zoom particles into the
            // output buffer

            // First, delete the low-res counterparts
            deleteParticles({Pos1,Pos2,Pos3,Vel1,Vel2,Vel3});


            // Then, inject the high-res particles
            insertParticles({Pos1,Pos2,Pos3,Vel1,Vel2,Vel3},{Pos1z,Pos2z,Pos3z,Vel1z,Vel2z,Vel3z});


            Mass = (MyFloat*)calloc(nPartTotal,sizeof(MyFloat));

            for(long i=0; i<nPartTotal; i++) {
                if(i<nPartLevel[0]-zoomParticleArray.size())
                    Mass[i] = pmass1;
                else
                    Mass[i] = pmass2;
            }

            if (gadgetformat==2){
                SaveGadget2( (base+ "zoom_gadget2.dat").c_str(), nPartTotal,
                CreateGadget2Header<MyFloat>(nPartTotal, 0, a, zin, boxlen[0], Om0, Ol0, hubble),
                Pos1, Vel1, Pos2, Vel2, Pos3, Vel3, Mass);
            }
            else if (gadgetformat==3){
                SaveGadget3( (base+ "zoom_gadget3.dat").c_str(), nPartTotal,
                CreateGadget3Header<MyFloat>(nPartTotal, 0, a, zin, boxlen[0], Om0, Ol0, hubble),
                Pos1, Vel1, Pos2, Vel2, Pos3, Vel3, Mass);
            }
            else if (gadgetformat==4) {
                SaveTipsy( (base+"zoom.tipsy").c_str(), nPartTotal, Pos1, Vel1, Pos2, Vel2, Pos3, Vel3,
			   boxlen[0],  Om0,  Ol0,  hubble,  a, pmass1, Mass, Ob0, 1.0/(140*n[0]));
            }


            free(Pos1);
            free(Vel1);
            free(Pos2);
            free(Vel2);
            free(Pos3);
            free(Vel3);

            dumpGrid(0);
            dumpGrid(1);

        }

    }

    virtual void prepare() {
        if(prepared)
            throw(std::runtime_error("Called prepare, but field is already prepared"));

        cerr << "***** PREPARE FIELD *****" << endl;
        prepared=true;
        UnderlyingField<MyFloat> *pField;

        for_each_level(level) nPartLevel[level] = ((long)n[level]*n[level])*n[level];

        base = make_base(indir);

        readCamb();
        drawRandom();

        if(n[1]>0)
            splitLevel0();

        applyPowerSpec();

        pField = new UnderlyingField<MyFloat>(this->P[0], this->pField_k[0], this->nPartLevel[0]);

        for_each_level(level) {
            pGrid[level] = new Grid<MyFloat>(n[level], dx[level], x_off[level], y_off[level], z_off[level]);
            if(level>0)
                pField->add_component(this->P[level],this->pField_k[level],this->nPartLevel[level]);
        }


        // TODO: actually combine the different levels before passing them to the constrainer
        pConstrainer = new MultiConstrainedField<MyFloat>(pField);
    }

    ///////////////////////
    // CONSTRAING CODE //
    ///////////////////////



protected:

    int deepestLevelWithParticles() {
        if(levelParticleArray[1].size()>0) return 1; else return 0;
    }

    MyFloat get_wrapped_delta(MyFloat x0, MyFloat x1) {
        MyFloat result = x0-x1;
        if(result>this->boxlen[0]/2) {
            result-=this->boxlen[0];
        }
        if(result<-this->boxlen[0]/2) {
            result+=this->boxlen[0];
        }
        return result;
    }


    void getCentre() {
        x0 = 0; y0 = 0; z0 =0;

        int level = deepestLevelWithParticles();

        MyFloat xa,ya,za, xb, yb, zb;
        pGrid[level]->get_centroid_location(levelParticleArray[level][0],xa,ya,za);

        for(long i=0;i<levelParticleArray[level].size();i++) {
            pGrid[level]->get_centroid_location(levelParticleArray[level][i],xb,yb,zb);
            x0+=get_wrapped_delta(xa,xb);
            y0+=get_wrapped_delta(ya,yb);
            z0+=get_wrapped_delta(za,zb);
        }
        x0/=levelParticleArray[level].size();
        y0/=levelParticleArray[level].size();
        z0/=levelParticleArray[level].size();
        x0+=xa;
        y0+=ya;
        z0+=za;
        cerr << "Centre of region is " << x0 << " " << y0 << " " << z0 << endl;
    }




    void AllocAndGetBuffer_int(const char* IDfile, bool append=false) {
        // count lines, allocate space, then read

        cerr << "Loading " << IDfile << endl;

        ifstream inf;
        inf.open(IDfile);

        if(!inf.is_open()) {
            throw(runtime_error("Error: could not open particle file"));
        }

        if(!append)
            genericParticleArray.clear();
        else
            cerr << "Append: starts with " << genericParticleArray.size() << " existing particles" << endl;


        while(true) {
            long x;
            inf >> x;
            if(inf.eof()) break;
            genericParticleArray.push_back(x);
        }

        cerr << "Total number of particles is " << genericParticleArray.size() << endl;

        interpretParticleList();
    }


    void cen_deriv4_alpha(long index, int direc, complex<MyFloat> *alpha, int level)

    {//4th order central difference

      MyFloat xp, yp, zp;//, zm2, zm1, zp2, zp1;
      pGrid[level]->get_centroid_location(index, xp, yp, zp);
      xp=get_wrapped_delta(xp,x0);
      yp=get_wrapped_delta(yp,y0);
      zp=get_wrapped_delta(zp,z0);

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

          // N.B. can't wrap - might be on subgrid

          ind_m1=pGrid[level]->find_next_ind_no_wrap(index, neg_step1);
          ind_p1=pGrid[level]->find_next_ind_no_wrap(index, step1);
          ind_m2=pGrid[level]->find_next_ind_no_wrap(ind_m1, neg_step1);
          ind_p2=pGrid[level]->find_next_ind_no_wrap(ind_p1, step1);

          MyFloat a=-1./12./dx[0], b=2./3./dx[0];  //the signs here so that L ~ - Nabla Phi

          alpha[ind_m2]+=(c[di]*a);
          alpha[ind_m1]+=(c[di]*b);
          alpha[ind_p1]+=(-c[di]*b);
          alpha[ind_p2]+=(-c[di]*a);
       }

    }

    complex<MyFloat> *calcConstraintVectorAllLevels(string name_in) {
        long nall = 0;
        for_each_level(level) nall+=nPartLevel[level];
        complex<MyFloat> *rval_k=(complex<MyFloat>*)calloc(nall,sizeof(complex<MyFloat>));
        for_each_level(level) calcConstraintVector(name_in, level, rval_k);
        return rval_k;
    }

    complex<MyFloat> *calcConstraintVector(string name_in, int level) {
        complex<MyFloat> *rval_k=(complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
        calcConstraintVector(name_in, level, rval_k, false);
        return rval_k;

    }

    void calcConstraintVector(string name_in, int level, complex<MyFloat> *ar, bool offset=true) {
        //
        // If offset is true, offset the storage to take account relative to the
        // start of ar, to take account of the other levels
        const char* name = name_in.c_str();
        complex<MyFloat> *rval=(complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
        complex<MyFloat> *rval_k;

        if(offset) {
            long offset_amount = 0;
            for(int l=0; l<level; l++) {
                offset_amount+=nPartLevel[l];
            }
            cerr << "Note level = " << level << " offset = " << offset_amount << endl;
            rval_k = &ar[offset_amount];
        } else {
            rval_k = ar;
        }

        if(strcasecmp(name,"overdensity")==0) {
            MyFloat w = 1.0/levelParticleArray[level].size();
            for(long i=0;i<levelParticleArray[level].size();i++) {
                rval[levelParticleArray[level][i]]+=w;
            }

            fft_r(rval_k, rval, this->n[level], 1);
        }
        else if(strcasecmp(name,"phi")==0) {
            MyFloat w = 1.0/levelParticleArray[level].size();
            for(long i=0;i<levelParticleArray[level].size();i++) {
                rval[levelParticleArray[level][i]]+=w;
            }
            complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
            fft_r(rval_kX, rval, this->n[level], 1);
            poiss(rval_k, rval_kX, this->n[level], this->boxlen[level], this->a, this->Om0);
            free(rval_kX);
        }
        else if(name[0]=='L' || name[0]=='l') {
            // angular momentum
            int direction=atoi(&name[1]);

            cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;

            for(long i=0;i<levelParticleArray[level].size();i++) {
                cen_deriv4_alpha(levelParticleArray[level][i], direction, rval, level);
            }
            complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
            fft_r(rval_kX, rval, this->n[level], 1);
            // The constraint as derived is on the potential. By considering
            // unitarity of FT, we can FT the constraint to get the constraint
            // on the density.
            poiss(rval_k, rval_kX, this->n[level], this->boxlen[level], this->a, this->Om0);
            free(rval_kX);
        } else {
            cout << "  -> UNKNOWN constraint vector type, returning zeros [this is bad]" << endl;

        }

        if(pField_x[level]!=NULL) {
            cout << " dot in real space = " << std::real(dot(pField_x[level], rval, this->nPartLevel[level])) << endl;
        }

        free(rval);
    }

    void freeRealDelta(int level=0) {
        if(pField_x[level]!=NULL) {
            free(pField_x[level]);
            pField_x[level]=NULL;
        }
    }

    void freeFourierDelta(int level=0) {
        if(pField_k[level]!=NULL) {
            free(pField_k[level]);
            pField_k[level]=NULL;
        }
    }

    void ensureRealDelta(int level=0) {
        if(pField_x[level]==NULL) {
            pField_x[level] = (complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
            fft_r(pField_x[level],this->pField_k[level],this->n[level],-1);
        }
    }

    void ensureFourierDelta(int level=0) {
        if(pField_k[level]==NULL) {
            pField_k[level] = (complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
            fft_r(pField_k[level],pField_x[level],n[level],1);
        }
    }

    void operateInRealSpace(int level=0) {
        ensureRealDelta(level);
        freeFourierDelta(level);
    }

    void operateInFourierSpace(int level=0) {
        ensureFourierDelta(level);
        freeRealDelta(level);
    }



public:


    void loadID(string fname) {
        AllocAndGetBuffer_int(fname.c_str());
        getCentre();
    }

    void appendID(string fname) {
        AllocAndGetBuffer_int(fname.c_str(), true);
        getCentre();
    }

    void centreParticle(long id) {
        pGrid[0]->get_centroid_location(id,x0,y0,z0);
    }

    void selectNearest() {
        MyFloat delta_x, delta_y, delta_z, r2_i;
        MyFloat xp,yp,zp;

        genericParticleArray.clear();

        for_each_level(level) {
            int repeat=1;
            if(level==0 && zoomfac>1) {
                repeat = zoomfac*zoomfac*zoomfac;
            }

            if(level==1) continue;

            levelParticleArray[level].clear();

            MyFloat r2_nearest = 1.0/0.0;
            long i_nearest;

            for(long i=0;i<this->nPartLevel[level];i++) {
                pGrid[level]->get_centroid_location(i,xp,yp,zp);
                delta_x = get_wrapped_delta(xp,x0);
                delta_y = get_wrapped_delta(yp,y0);
                delta_z = get_wrapped_delta(zp,z0);
                r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;

                if(r2_i<r2_nearest) {
                    r2_nearest = r2_i;
                    i_nearest = i;
                }

            }

            for(int q=0; q<repeat; q++)
                levelParticleArray[level].push_back(i_nearest);

        }

    }

    void selectSphere(float radius) {
        MyFloat r2 = radius*radius;
        MyFloat delta_x, delta_y, delta_z, r2_i;
        MyFloat xp,yp,zp;

        genericParticleArray.clear();

        for_each_level(level) {
            int repeat=1;
            if(level==0 && zoomfac>1) {
                repeat = zoomfac*zoomfac*zoomfac;
            }


            levelParticleArray[level].clear();
            for(long i=0;i<this->nPartLevel[level];i++) {
                pGrid[level]->get_centroid_location(i,xp,yp,zp);
                delta_x = get_wrapped_delta(xp,x0);
                delta_y = get_wrapped_delta(yp,y0);
                delta_z = get_wrapped_delta(zp,z0);
                r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
                if(r2_i<r2)
                    for(int q=0; q<repeat; q++)
                        levelParticleArray[level].push_back(i);
            }
        }

    }


    void centreDenmax() {

        float den_max=-1000;
        long index_max=0;
        ensureRealDelta();

        for(long i=0;i<genericParticleArray.size();i++) {
            if(std::real(pField_x[0][genericParticleArray[i]])>den_max) {
                index_max=genericParticleArray[i];
                den_max = std::real(pField_x[0][genericParticleArray[i]]);
                pGrid[0]->get_centroid_location(i,x0,y0,z0);
            }
        }

        cerr << "Denmax = " << den_max <<", index=" << index_max << " coords=" << x0 << " " << y0 << " " << z0 << endl;

    }

    void setCentre(MyFloat xin, MyFloat yin, MyFloat zin) {
            x0=xin;
            y0=yin;
            z0=zin;
    }

    void reorderBuffer() {

        cout << "Reordering buffer radially..." << endl;
        cout << " [taking centre = " << x0 << " " << y0 << " " <<z0 << "]" << endl;

        std::vector<MyFloat> r2(genericParticleArray.size());
        MyFloat delta_x, delta_y, delta_z;

        std::vector<size_t> index(genericParticleArray.size());

        for(int i=0;i<genericParticleArray.size();i++) {
            pGrid[0]->get_centroid_location(i,delta_x, delta_y, delta_z);
            delta_x = get_wrapped_delta(delta_x,x0);
            delta_y = get_wrapped_delta(delta_y,y0);
            delta_z = get_wrapped_delta(delta_z,z0);
            r2[i] = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
            index[i]=i;
        }

        // Now sort the index array
        std::sort(index.begin(),index.end(),
             [&r2](size_t i1, size_t i2) { return r2[i1]<r2[i2]; } );

        // Turn the index array into something pointing to the particles
        for(int i=0; i<genericParticleArray.size(); i++) {
            index[i] = genericParticleArray[index[i]];
        }

        // Copy back into the particle array
        for(int i=0; i<genericParticleArray.size(); i++) {
            genericParticleArray[i] = index[i];
        }

    }

    void truncateBuffer(float x) {
        if(x<0 ) {
            cerr << "Truncate command takes a fraction between 0 and 1 or a number>=1" << endl;
            exit(0);
        }
        if(x<1)
            genericParticleArray.resize(((int)(genericParticleArray.size()*x)));
        else
            genericParticleArray.resize(((int)x));
    }

    void calculate(string name) {
        float val=0.0;
        for_each_level(level) {
            complex<MyFloat> *vec = calcConstraintVector(name,level);
            val+=std::real(dot(vec, this->pField_k[level], this->nPartLevel[level]));
            cerr << val << " ";
            free(vec);
        }
        cerr << endl;
        cout << name << ": calculated value = " <<  val << endl;
    }

    void constrain(string name, string type, float value) {
        if(pConstrainer==NULL) {
            throw runtime_error("No constraint information is available. Is your done command too early, or repeated?");
        }

        bool relative=false;
        if (strcasecmp(type.c_str(),"relative")==0) {
            relative=true;
        } else if (strcasecmp(type.c_str(),"absolute")!=0) {
            cerr << "Constraints must state either relative or absolute" << endl;
            exit(0);
        }

        std::complex<MyFloat> constraint = value;
        complex<MyFloat> *vec = calcConstraintVectorAllLevels(name);
        std::complex<MyFloat> initv = pConstrainer->underlying->v1_dot_y(vec);

        if(relative) constraint*=initv;

        cout << name << ": initial value = " << initv << ", constraining to " << constraint << endl;
        pConstrainer->add_constraint(vec, constraint, initv);

    }

    void cov() {
        if(pConstrainer==NULL) {
            throw runtime_error("No constraint information is available. Is your done command too early, or repeated?");
        }
        pConstrainer->print_covariance();
    }


    void done() {
        if(pConstrainer==NULL) {
            throw runtime_error("No constraint information is available. Is your done command too early, or repeated?");
        }

        long nall = 0;
        for_each_level(level) nall+=nPartLevel[level];
        complex<MyFloat> *y1=(complex<MyFloat>*)calloc(nall,sizeof(complex<MyFloat>));

        pConstrainer->prepare();
        pConstrainer->get_realization(y1);

        long j=0;
        for_each_level(level) {
            for(long i=0; i<this->nPartLevel[level]; i++) {
                this->pField_k[level][i]=y1[j];
                j++;
            }
        }

        free(y1);

        cout << "Expected Delta chi^2=" << pConstrainer->get_delta_chi2() << endl;

        delete pConstrainer;
        pConstrainer=NULL;


        /*
        // the following wants to come JUST before writing...
        // and perhaps should be performed after k-filtering to get the
        // directional offsets (that's where it can currently be found)
        for_each_level(level) {
            if(level>0)
                interpolateLevel(level); // copy in the underlying field
        }

        recombineLevel0();
        */


    }

    void reverse() {
        if(pConstrainer!=NULL) {
            throw runtime_error("Can only reverse field direction after a 'done' command finailises the constraints");
        }

        for(long i=0; i<this->nPartLevel[0]; i++)
            this->pField_k[0][i]=-this->pField_k[0][i];
    }


};
