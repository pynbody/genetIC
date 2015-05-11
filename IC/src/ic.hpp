#include <string>
#include <tuple>
#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <cassert>
#include <functional>
#include <algorithm>
#include <memory>

#include "numpy.hpp"

#define for_each_level(level) for(int level=0; level<2 && n[level]>0; level++)

using namespace std;


template<typename MyFloat>
class IC {
protected:



    MyFloat Om0, Ol0, Ob0, hubble, zin, a, sigma8, ns;

    // Everything about the grids:
    MyFloat boxlen[2], dx[2];      // the box length of each grid
    int n[2];                      // number of grid divisions along one axis

    std::vector<std::shared_ptr<Grid<MyFloat>>> pGrid;       // the objects that help us relate points on the grid



    complex<MyFloat> *P[2]; // the power spectrum for each k-space point

    std::vector<complex<MyFloat>> pField_k_0_high; // high-k modes for level 0

    MyFloat x_off[2], y_off[2], z_off[2]; // x,y,z offsets for subgrids

    int zoomfac; // the zoom factor, i.e. the ratio dx/dx[1]
    MyFloat k_cut;   // the cut-off k for low-scale power


    int out, gadgetformat, seed;


    int nCambLines;

    CAMB<MyFloat> spectrum;

    size_t nPartLevel[2];

    string incamb, indir, inname, base;

    bool whiteNoiseFourier;

    bool prepared;




    std::vector<size_t> genericParticleArray;



    MyFloat x0, y0, z0;

    MultiConstrainedField<MyFloat> *pConstrainer;
    shared_ptr<ParticleMapper<MyFloat>> pMapper;

    using RefFieldType = decltype(pGrid[0]->get_field());


public:
    IC() : pMapper(new ParticleMapper<MyFloat>())
    {
        pConstrainer=NULL;
        whiteNoiseFourier=false;
        hubble=0.701;   // old default
        Ob0=-1.0;
        ns = 0.96;      // old default
        n[0]=-1;
        n[1]=-1; // no subgrid by default
        boxlen[0]=-1;
        x_off[0]=y_off[0]=z_off[0]=0;
        x_off[1]=y_off[1]=z_off[1]=0;
        prepared = false;
    }

    ~IC() {

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
        initGrid();
    }

    void setZ0(MyFloat in) {
        zin = in;
        a=1./(zin+1.);
    }

    void setn(int in) {
        n[0] = in;
        initGrid();
    }

    void initGrid(int level=0) {
        if(n[level]<0 || boxlen[level]<0)
            return;

        nPartLevel[level] = ((long)n[level]*n[level])*n[level];
        dx[level] = boxlen[level]/n[level];

        if(pGrid.size()!=level)
            throw std::runtime_error("Trying to re-initialize a grid level");

        pGrid.emplace_back(new Grid<MyFloat>(boxlen[0], n[level], dx[level],
                                             x_off[level], y_off[level], z_off[level]));

        if(level==0)
            pMapper = std::shared_ptr<ParticleMapper<MyFloat>>(new OneLevelParticleMapper<MyFloat>(pGrid.back()));
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

        AllocAndGetBuffer_int(fname);
        auto zoomParticleArray = pGrid[0]->particleArray;

        // Sorting now happens inside mapper class
        // std::sort(zoomParticleArray.begin(), zoomParticleArray.end());



        // find boundaries
        int x0, x1, y0, y1, z0, z1;
        int x,y,z;

        x0=y0=z0=n[0];
        x1=y1=z1=0;
        Grid<MyFloat> g(n[0]);

        // TO DO: wrap the box sensibly

        for(size_t i=0; i<zoomParticleArray.size(); i++) {
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



        initZoom(zoomParticleArray);

    }

    void initZoom(const std::vector<size_t> & zoomParticleArray) {
        zoomfac = (boxlen[0]/n[0])/(boxlen[1]/n[1]);
        dx[1] = boxlen[1]/n[1];
        cout << "Initialized a zoom region:" << endl;
        cout << "  Subbox length   = " << boxlen[1] << " Mpc/h" << endl;
        cout << "  n[1]            = " << n[1] << endl;
        cout << "  dx[1]           = " << dx[1] << endl;
        cout << "  Zoom factor     = " << zoomfac << endl;
        cout << "  Low-left corner = " << x_off[1] << ", " << y_off[1] << ", " << z_off[1] << endl;

        initGrid(1);

        auto pMapper1 = std::shared_ptr<ParticleMapper<MyFloat>>(new OneLevelParticleMapper<MyFloat>(pGrid.back()));

        pMapper = std::shared_ptr<ParticleMapper<MyFloat>>(
            new TwoLevelParticleMapper<MyFloat>(pMapper, pMapper1, zoomParticleArray, zoomfac*zoomfac*zoomfac));


        cout << "  Total particles = " << pMapper->size() << endl;
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
        spectrum.read(incamb);
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
                                  MyFloat norm, int n, RefFieldType pField_k) {
        long id_k, id_negk;
        id_k = kgridToIndex(k1,k2,k3,n);
        id_negk = kgridToIndex(-k1,-k2,-k3,n);
        pField_k[id_k]=std::complex<MyFloat>(norm*gsl_ran_gaussian_ziggurat(r,1.),norm*gsl_ran_gaussian_ziggurat(r,1.));

        // reality condition:
        pField_k[id_negk]=std::conj(pField_k[id_k]);
    }


    static void drawRandomForSpecifiedGrid(int n,
                                           RefFieldType pField_k,
                                           gsl_rng * r)
    {
        long nPartTotal = n;
        nPartTotal*=n*n;

        MyFloat sigma=sqrt((MyFloat)(nPartTotal));

        complex<MyFloat> *rnd=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

        cerr << "Drawing random numbers..."<< endl;

        long i;

        // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
        for(i=0;i<nPartTotal;i++){rnd[i]=gsl_ran_gaussian_ziggurat(r,1.)*sigma;}// cout<< "rnd "<< rnd[i] << endl;}


        cout<< "First FFT..." <<endl;
        fft(pField_k.data(), rnd, n, 1);

        free(rnd);

        pField_k[0]=complex<MyFloat>(0.,0.); //assure mean==0



    }

    static void drawRandomForSpecifiedGridFourier(int n,
                                           RefFieldType pField_k,
                                           gsl_rng * r) {

       long nPartTotal = n;
       nPartTotal*=n*n;

       MyFloat sigma=sqrt((MyFloat)(nPartTotal));


       // do it in fourier space, in order of increasing |k|, so that
       // resolution can be scaled by factors of 2 and we still get
       // the 'same' field

       cerr << "Drawing random numbers in fourier space..."<< endl;
       int ks,k1,k2;

       sigma/=sqrt(2.0);

       // N.B. DO NOT PARALLELIZE this loop - want things to be done in a reliable order
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
                drawRandomForSpecifiedGridFourier(n[level], pGrid[level]->get_field_fourier(),r);
            else
                drawRandomForSpecifiedGrid(n[level], pGrid[level]->get_field_fourier(),r);
        }

        gsl_rng_free (r);
    }

    void zeroLevel(int level) {
        cerr << "*** Warning: your script calls zeroLevel("<<level <<"). This is intended for testing purposes only!" << endl;
        if(level==-1) {
            std::fill(pField_k_0_high.begin(), pField_k_0_high.end(), 0);
        } else {
            auto &field = pGrid[level]->get_field_real();
            std::fill(field.begin(), field.end(), 0);
        }
    }

    template<typename T>
    void interpolateLevel(int level, T* pField_x_l, T* pField_x_p) {

        if(level<=0)
            throw std::runtime_error("Trying to interpolate onto the top-level grid");

        int n_l = n[level];
        MyFloat dx_l = dx[level];
        MyFloat dx_p = dx[level-1];

        auto pGrid_l = pGrid[level];
        auto pGrid_p = pGrid[level-1];


        // get offsets of bottom-left of grids
        MyFloat x_off_l = x_off[level], y_off_l = y_off[level], z_off_l = z_off[level];
        MyFloat x_off_p = x_off[level-1], y_off_p = y_off[level-1], z_off_p = z_off[level-1];


        #pragma omp parallel for schedule(static)
        for(int x_l=0; x_l<n_l; x_l++) {

            // private variables used in calculation:

            MyFloat x,y,z; // coordinates
            int x_p_0, y_p_0, z_p_0, x_p_1, y_p_1, z_p_1; // parent grid locations
            MyFloat xw0,yw0,zw0; // weights for lower-left parent grid contribution
            MyFloat xw1,yw1,zw1; // weights for upper-right parent grid contribution


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
        RefFieldType pField_x_l = pGrid[level]->get_field_real();
        RefFieldType pField_x_p = pGrid[level-1]->get_field_real();

        interpolateLevel(level,pField_x_l.data(),pField_x_p.data());

        auto offset_fields_l = pGrid[level]->get_offset_fields();
        auto offset_fields_p = pGrid[level-1]->get_offset_fields();

        if(std::get<0>(offset_fields_l)->size()>0 && std::get<0>(offset_fields_p)->size()>0) {
            interpolateLevel(level,std::get<0>(offset_fields_l)->data(),std::get<0>(offset_fields_p)->data());
            interpolateLevel(level,std::get<1>(offset_fields_l)->data(),std::get<1>(offset_fields_p)->data());
            interpolateLevel(level,std::get<2>(offset_fields_l)->data(),std::get<2>(offset_fields_p)->data());
        }
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
        pField_k_0_high = pGrid[0]->get_field_fourier();
        cerr << "Apply transfer function for high-k part of level 0 field" << endl;
        cerr << "  sizes = " << pGrid[0]->get_field_fourier().size() << " " <<pField_k_0_high.size() << endl;
        cerr << "  first el = " << pField_k_0_high[10] << " " << pGrid[0]->get_field_fourier()[10] << endl;
        applyPowerSpecForLevel(0,true);
        cerr << "  first el = " << pField_k_0_high[10] << " " << pGrid[0]->get_field_fourier()[10] << endl;
    }

    void recombineLevel0() {

        RefFieldType pField_k = pGrid[0]->get_field_fourier();

        for(size_t i=0; i<this->nPartLevel[0]; i++)
            pField_k[i]+=pField_k_0_high[i];

        pField_k_0_high.clear();

    }

    void applyPowerSpecForLevel(int level, bool high_k=false) {
        //scale white-noise delta with initial PS

        MyFloat grwfac;

        //growth factor normalised to 1 today:
        grwfac=D(a, Om0, Ol0)/D(1., Om0, Ol0);

        cout<< "Growth factor " << grwfac << endl;
        MyFloat sg8;

        sg8=spectrum.sig(8., ns, boxlen[level], n[level]);
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

        // TODO - fix up types here

        auto & pField_k_this = high_k?pField_k_0_high:pGrid[level]->get_field_fourier();


        spectrum.applyTransfer(n[level], kw, ns, norm_amp,
                               pField_k_this, pField_k_this, P[level], filter);

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

        ostringstream filename;
        filename << indir << "/grid-" << level << ".npy";

        const int dim[3] = { n[level],n[level],n[level] };
        numpy::SaveArrayAsNumpy(filename.str(), true, 3, dim, pGrid[level]->get_field_real().data() );


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
        powsp_noJing(n[level], pGrid[level]->get_field_real().data(),
                     (base+"_"+((char)(level+'0'))+".ps").c_str(), boxlen[level]);
    }

    void zeldovich(int level) {
        //Zeldovich approx.

        MyFloat hfac=1.*100.*sqrt(Om0/a/a/a+Ol0)*sqrt(a);
        //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
        //TODO: hardcoded value of f=1 is inaccurate, but fomega currently gives wrong nults

        MyFloat pmass=27.78*Om0*powf(boxlen[level]/(MyFloat)(n[level]),3.0);

        pGrid[level]->zeldovich(hfac,pmass);

    }



    void zeldovich() {
        if(n[1]<=0) {
            cerr << "***** WRITE - no zoom *****" << endl;
            zeldovich(0);
        } else {
            zeldovich(0);
            zeldovich(1);

            cerr << "Interpolating low-frequency information into zoom region...";
            cerr.flush();
            interpolateLevel(1);
            cerr << "done." << endl;


            cerr << "Re-introducing high-k modes into low-res region...";
            recombineLevel0();
            zeldovich(0);
            cerr << "done." << endl;

        }
    }

    ///////////////////////////////////////////////
    // ZOOM particle management
    ///////////////////////////////////////////////



    void interpretParticleList() {
        // copies the overall particle list into the relevant grid levels
        pMapper->interpretParticleList(genericParticleArray);
    }



    void write() {

        zeldovich();

        auto finalMapper=pMapper;

        if(Ob0>0) {


            // Add gas only to the deepest level. Pass the whole pGrid
            // vector if you want to add gas to every level.
            auto gasMapper = pMapper->addGas(Ob0/Om0,
                                            {pGrid.back()});


            // graft the gas particles onto the start of the map
            finalMapper = std::make_shared<AddGasMapper<MyFloat>>(
                gasMapper, pMapper, true);


        }

        cerr << "Write, ndm=" << finalMapper->size_dm() << ", ngas=" << finalMapper->size_gas() << endl;

        /*
        if (gadgetformat==2){
            SaveGadget2( (base+ ".gadget2").c_str(), nPartTotal,
            CreateGadget2Header<MyFloat>(nPartTotal, 0, a, zin, boxlen[0], Om0, Ol0, hubble),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3, Mass);
        }
        else if (gadgetformat==3){
            SaveGadget3( (base+ ".gadget3").c_str(), nPartTotal,
            CreateGadget3Header<MyFloat>(nPartTotal, 0, a, zin, boxlen[0], Om0, Ol0, hubble),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3, Mass);
        }
        */
        if (gadgetformat==4)
            SaveTipsy(base+".tipsy", boxlen[0], Om0, Ol0, hubble, a, finalMapper);



    }

    virtual void prepare() {
        if(prepared)
            throw(std::runtime_error("Called prepare, but grid is already prepared for constraints"));

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

        pField = new UnderlyingField<MyFloat>(this->P[0], this->pGrid[0], this->nPartLevel[0]);

        for_each_level(level) {
            if(level>0)
                pField->add_component(this->P[level],
                                      this->pGrid[level],
                                      this->nPartLevel[level]);
        }

        // TODO: actually combine the different levels before passing them to the constrainer
        pConstrainer = new MultiConstrainedField<MyFloat>(pField);
    }

    ///////////////////////
    // CONSTRAING CODE //
    ///////////////////////



protected:

    int deepestLevelWithParticles() {
        if(pGrid.size()>1 && pGrid[1]->particleArray.size()>0) return 1; else return 0;
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
        pGrid[level]->get_centroid_location(pGrid[level]->particleArray[0],xa,ya,za);

        for(size_t i=0;i<pGrid[level]->particleArray.size();i++) {
            pGrid[level]->get_centroid_location(pGrid[level]->particleArray[i],xb,yb,zb);
            x0+=get_wrapped_delta(xa,xb);
            y0+=get_wrapped_delta(ya,yb);
            z0+=get_wrapped_delta(za,zb);
        }
        x0/=pGrid[level]->particleArray.size();
        y0/=pGrid[level]->particleArray.size();
        z0/=pGrid[level]->particleArray.size();
        x0+=xa;
        y0+=ya;
        z0+=za;
        cerr << "Centre of region is " << x0 << " " << y0 << " " << z0 << endl;
    }




    void AllocAndGetBuffer_int(std::string IDfile, bool append=false) {

        cerr << "Loading " << IDfile << endl;

        if(!append)
            genericParticleArray.clear();
        else
            cerr << " -> append: starts with " << genericParticleArray.size() << " existing particles" << endl;

        getBuffer(genericParticleArray, IDfile);

        cerr << "  -> total number of particles is " << genericParticleArray.size() << " " << genericParticleArray[0] << " " << genericParticleArray.back() << endl;

        interpretParticleList();
    }


    void cen_deriv4_alpha(long index, int direc, complex<MyFloat> *alpha, int level)

    {//4th order central difference

      MyFloat xp, yp, zp;//, zm2, zm1, zp2, zp1;
      pGrid[level]->get_centroid_location(index, xp, yp, zp);
      xp=get_wrapped_delta(xp,x0);
      yp=get_wrapped_delta(yp,y0);
      zp=get_wrapped_delta(zp,z0);


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
            MyFloat w = 1.0/pGrid[level]->particleArray.size();
            for(size_t i=0;i<pGrid[level]->particleArray.size();i++) {
                rval[pGrid[level]->particleArray[i]]+=w;
            }

            fft(rval_k, rval, this->n[level], 1);
        }
        else if(strcasecmp(name,"phi")==0) {
            MyFloat w = 1.0/pGrid[level]->particleArray.size();
            for(size_t i=0;i<pGrid[level]->particleArray.size();i++) {
                rval[pGrid[level]->particleArray[i]]+=w;
            }
            complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
            fft(rval_kX, rval, this->n[level], 1);
            poiss(rval_k, rval_kX, this->n[level], this->boxlen[level], this->a, this->Om0);
            free(rval_kX);
        }
        else if(name[0]=='L' || name[0]=='l') {
            // angular momentum
            int direction=atoi(&name[1]);

            cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;

            for(size_t i=0;i<pGrid[level]->particleArray.size();i++) {
                cen_deriv4_alpha(pGrid[level]->particleArray[i], direction, rval, level);
            }
            complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(this->nPartLevel[level],sizeof(complex<MyFloat>));
            fft(rval_kX, rval, this->n[level], 1);
            // The constraint as derived is on the potential. By considering
            // unitarity of FT, we can FT the constraint to get the constraint
            // on the density.
            poiss(rval_k, rval_kX, this->n[level], this->boxlen[level], this->a, this->Om0);
            free(rval_kX);
        } else {
            cout << "  -> UNKNOWN constraint vector type, returning zeros [this is bad]" << endl;

        }


        // cout << " dot in real space = " << std::real(dot(pField_x[level], rval, this->nPartLevel[level])) << endl;


        free(rval);
    }




public:


    void loadID(string fname) {
        AllocAndGetBuffer_int(fname);
        getCentre();
    }

    void appendID(string fname) {
        AllocAndGetBuffer_int(fname, true);
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

            pGrid[level]->particleArray.clear();

            MyFloat r2_nearest = 1.0/0.0;
            long i_nearest;

            for(size_t i=0;i<this->nPartLevel[level];i++) {
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
                pGrid[level]->particleArray.push_back(i_nearest);

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


            pGrid[level]->particleArray.clear();
            for(size_t i=0;i<this->nPartLevel[level];i++) {
                pGrid[level]->get_centroid_location(i,xp,yp,zp);
                delta_x = get_wrapped_delta(xp,x0);
                delta_y = get_wrapped_delta(yp,y0);
                delta_z = get_wrapped_delta(zp,z0);
                r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
                if(r2_i<r2)
                    for(int q=0; q<repeat; q++)
                        pGrid[level]->particleArray.push_back(i);
            }
        }

    }


    void centreDenmax() {

        float den_max=-1000;
        long index_max=0;

        auto & pField_x = pGrid[0]->get_field_real();
        for(size_t i=0;i<genericParticleArray.size();i++) {
            if(std::real(pField_x[genericParticleArray[i]])>den_max) {
                index_max=genericParticleArray[i];
                den_max = std::real(pField_x[genericParticleArray[i]]);
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

        for(size_t i=0;i<genericParticleArray.size();i++) {
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
        for(size_t i=0; i<genericParticleArray.size(); i++) {
            index[i] = genericParticleArray[index[i]];
        }

        // Copy back into the particle array
        for(size_t i=0; i<genericParticleArray.size(); i++) {
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
            val+=std::real(dot(vec, pGrid[level]->get_field_fourier().data(), this->nPartLevel[level]));
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
            for(size_t i=0; i<this->nPartLevel[level]; i++) {
                this->pGrid[level]->get_field_fourier()[i]=y1[j];
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


        // All done - write out
        write();



    }

    void reverse() {
        if(pConstrainer!=NULL) {
            throw runtime_error("Can only reverse field direction after a 'done' command finailises the constraints");
        }

        for_each_level(level) {
            auto & field = pGrid[level]->get_field();
            for(size_t i=0; i<this->nPartLevel[0]; i++)
                field[i]=-field[i];
        }
    }


};
