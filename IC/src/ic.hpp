#ifndef _IC_HPP_INCLUDED
#define _IC_HPP_INCLUDED

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
#include <limits>

#include "numpy.hpp"

#define for_each_level(level) for(int level=0; level<2 && n[level]>0; level++)

using namespace std;

template<typename MyFloat>
class DummyIC;

template<typename MyFloat>
class IC {
protected:


    friend class DummyIC<MyFloat>;

    MyFloat Om0, Ol0, Ob0, hubble, zin, a, sigma8, ns;

    // Everything about the grids:
    MyFloat boxlen[2], dx[2];      // the box length of each grid
    int n[2];                      // number of grid divisions along one axis

    int supersample, subsample;               // DM supersampling to perform on zoom grid, and subsampling on base grid

    std::vector<std::shared_ptr<Grid<MyFloat>>> pGrid;       // the objects that help us relate points on the grid



    complex<MyFloat> *P[2]; // the power spectrum for each k-space point

    std::vector<complex<MyFloat>> pField_k_0_high; // high-k modes for level 0

    MyFloat x_off[2], y_off[2], z_off[2]; // x,y,z offsets for subgrids

    MyFloat xOffOutput, yOffOutput, zOffOutput;

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
    std::vector<size_t> zoomParticleArray;


    MyFloat x0, y0, z0;

    MultiConstrainedField<MyFloat> *pConstrainer;
    shared_ptr<ParticleMapper<MyFloat>> pMapper;
    shared_ptr<ParticleMapper<MyFloat>> pInputMapper;

    using RefFieldType = decltype(pGrid[0]->getField());

    ClassDispatch<IC<MyFloat>, void> &interpreter;


public:
    IC(ClassDispatch<IC<MyFloat>, void> &interpreter) : interpreter(interpreter), pMapper(new ParticleMapper<MyFloat>())
    {
        pConstrainer=nullptr;
        pInputMapper = nullptr;
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
        supersample = 1;
        subsample = 1;
        xOffOutput = 0;
        yOffOutput = 0;
        zOffOutput = 0;
    }

    ~IC() {

    }

    void setOmegaM0(MyFloat in) {
        Om0=in;
    }

    void setOmegaB0(MyFloat in) {
        Ob0=in;
        // now that we have gas, mapper may have changed:
        initMapper();
    }

    void setOmegaLambda0(MyFloat in) {
        Ol0=in;
    }

    void setHubble(MyFloat in) {
        hubble=in;
    }

    void offsetOutput(MyFloat x, MyFloat y, MyFloat z) {
        xOffOutput = x;
        yOffOutput = y;
        zOffOutput = z;
        initMapper();
    }

    void setSigma8(MyFloat in) {
        sigma8 = in;
    }

    void setSupersample(int in) {
        supersample = in;
        initMapper();
    }

    void setSubsample(int in) {
        subsample = in;
        initMapper();
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

    virtual void initGrid(unsigned int level=0) {

        if(n[level]<0 || boxlen[level]<0)
            return;

        nPartLevel[level] = ((long)n[level]*n[level])*n[level];
        dx[level] = boxlen[level]/n[level];

        if(pGrid.size()!=level)
            throw std::runtime_error("Trying to re-initialize a grid level");

        pGrid.push_back(std::make_shared<Grid<MyFloat>>(boxlen[0], n[level], dx[level],
                                             x_off[level], y_off[level], z_off[level]));

        // new level, particle mapper has changed:
        initMapper();

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
      AllocAndGetBuffer_int(fname);
      doZoom();
    }

    void doZoom() {
        if(n[1]==0)
            throw(std::runtime_error("Set n2 before specifying the zoom particles"));

        if(boxlen[1]==0)
            throw(std::runtime_error("Set the zoom factor before specifying the zoom particles"));

        pGrid[0]->gatherParticleList(zoomParticleArray);

        // find boundaries
        int x0, x1, y0, y1, z0, z1;
        int x,y,z;

        x0=y0=z0=n[0];
        x1=y1=z1=0;
        Grid<MyFloat> g(n[0]);

        // TO DO: wrap the box sensibly

        for(size_t i=0; i<zoomParticleArray.size(); i++) {
            g.getCoordinates(zoomParticleArray[i],x,y,z);
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
        cout << "Initialized a zoom region:" << endl;
        cout << "  Subbox length   = " << boxlen[1] << " Mpc/h" << endl;
        cout << "  n[1]            = " << n[1] << endl;
        cout << "  dx[1]           = " << dx[1] << endl;
        cout << "  Zoom factor     = " << zoomfac << endl;
        cout << "  Low-left corner = " << x_off[1] << ", " << y_off[1] << ", " << z_off[1] << endl;
        cout << "  Num particles   = " << zoomParticleArray.size() << endl;
        initGrid(1);

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

        // these need to be initialized in explicit order - can't leave it to compilers
        // to choose as they choose differently...
        MyFloat a = norm*gsl_ran_gaussian_ziggurat(r,1.);
        MyFloat b = norm*gsl_ran_gaussian_ziggurat(r,1.);

        pField_k[id_k]=std::complex<MyFloat>(a,b);

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

    virtual void drawRandom() {
        gsl_rng * r;
        const gsl_rng_type * T; //generator type variable

        T = gsl_rng_ranlxs2; // shouldn't this be gsl_rng_ranlxd2 for MyFloat = double?
        r = gsl_rng_alloc (T); //this allocates memory for the generator with type T
        gsl_rng_set(r,seed);

        for_each_level(level) {
            cerr << "Generate random noise for level " << level <<endl;
            if(whiteNoiseFourier)
                drawRandomForSpecifiedGridFourier(n[level], pGrid[level]->getFieldFourier(),r);
            else
                drawRandomForSpecifiedGrid(n[level], pGrid[level]->getFieldFourier(),r);
        }

        gsl_rng_free (r);
    }

    virtual void zeroLevel(int level) {
        cerr << "*** Warning: your script calls zeroLevel("<<level <<"). This is intended for testing purposes only!" << endl;
        if(level==-1) {
            std::fill(pField_k_0_high.begin(), pField_k_0_high.end(), 0);
        } else {
            auto &field = pGrid[level]->getFieldReal();
            std::fill(field.begin(), field.end(), 0);
        }
    }

    template<typename T>
    void interpolateIntoLevel(int level, T& pField_x_l, T& pField_x_p) {

        if(level<=0)
            throw std::runtime_error("Trying to interpolate onto the top-level grid");

        auto pGrid_l = pGrid[level];
        auto pGrid_p = pGrid[level-1];

        assert(pField_x_l.size()==pGrid_l->size3);
        assert(pField_x_p.size()==pGrid_p->size3);


        #pragma omp parallel for schedule(static)
        for(size_t ind_l=0; ind_l<pGrid_l->size3; ind_l++) {
            MyFloat x,y,z;
            pGrid_l->getCentroidLocation(ind_l,x,y,z);
            pField_x_l[ind_l]+=pGrid_p->getFieldInterpolated(x,y,z,pField_x_p);
        }

    }

    virtual void interpolateIntoLevel(int level) {
        RefFieldType pField_x_l = pGrid[level]->getFieldReal();
        RefFieldType pField_x_p = pGrid[level-1]->getFieldReal();

        interpolateIntoLevel(level,pField_x_l,pField_x_p);

        auto offset_fields_l = pGrid[level]->getOffsetFields();
        auto offset_fields_p = pGrid[level-1]->getOffsetFields();

        if(std::get<0>(offset_fields_l)->size()>0 && std::get<0>(offset_fields_p)->size()>0) {
            interpolateIntoLevel(level,*std::get<0>(offset_fields_l),*std::get<0>(offset_fields_p));
            interpolateIntoLevel(level,*std::get<1>(offset_fields_l),*std::get<1>(offset_fields_p));
            interpolateIntoLevel(level,*std::get<2>(offset_fields_l),*std::get<2>(offset_fields_p));
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


    virtual void splitLevel0() {
        pField_k_0_high = pGrid[0]->getFieldFourier();
        cerr << "Apply transfer function for high-k part of level 0 field" << endl;
        cerr << "  sizes = " << pGrid[0]->getFieldFourier().size() << " " <<pField_k_0_high.size() << endl;
        cerr << "  first el = " << pField_k_0_high[10] << " " << pGrid[0]->getFieldFourier()[10] << endl;
        applyPowerSpecForLevel(0,true);
        cerr << "  first el = " << pField_k_0_high[10] << " " << pGrid[0]->getFieldFourier()[10] << endl;
    }

    virtual void recombineLevel0() {

        RefFieldType pField_k = pGrid[0]->getFieldFourier();

        for(size_t i=0; i<this->nPartLevel[0]; i++)
            pField_k[i]+=pField_k_0_high[i];

        pField_k_0_high.clear();

    }

    virtual void applyPowerSpecForLevel(int level, bool high_k=false) {
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

        auto & pField_k_this = high_k?pField_k_0_high:pGrid[level]->getFieldFourier();


        spectrum.applyTransfer(n[level], kw, ns, norm_amp,
                               pField_k_this, pField_k_this, P[level], filter);

        //assert(abs(real(norm_iter)) >1e-12); //norm!=0
        //assert(abs(imag(norm_iter)) <1e-12); //because there shouldn't be an imaginary part since we assume d is purely real
        //TODO think about this (it fails more often than it should?)

        cout<<"Transfer applied!"<<endl;
        // cout<< "Power spectrum sample: " << P[level][0] << " " << P[level][1] <<" " << P[level][nPartTotal-1] <<endl;

        // std::cout<<"Initial chi^2 (white noise, fourier space) = " << chi2(pField_k[level],P[level],nPartTotal) << std::endl;



    }

    virtual void applyPowerSpec() {
        for_each_level(level) {
            cerr << "Apply transfer function for level " << level <<endl;
            applyPowerSpecForLevel(level);
        }
    }

    virtual void dumpGrid(int level=0) {

        ostringstream filename;
        filename << indir << "/grid-" << level << ".npy";

        const int dim[3] = { n[level],n[level],n[level] };
        numpy::SaveArrayAsNumpy(filename.str(), true, 3, dim, pGrid[level]->getFieldReal().data() );


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

    virtual void dumpPS(int level=0) {
        powsp_noJing(n[level], pGrid[level]->getFieldReal().data(),
                     (base+"_"+((char)(level+'0'))+".ps").c_str(), boxlen[level]);
    }

    virtual void zeldovichForLevel(int level) {
        //Zeldovich approx.

        MyFloat hfac=1.*100.*sqrt(Om0/a/a/a+Ol0)*sqrt(a);
        //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
        //TODO: hardcoded value of f=1 is inaccurate, but fomega currently gives wrong nults

        MyFloat pmass=27.78*Om0*powf(boxlen[level]/(MyFloat)(n[level]),3.0);

        pGrid[level]->zeldovich(hfac,pmass);

    }



    virtual void zeldovich() {
        if(pGrid.size()==0) {
          throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
        } else if(pGrid.size()==1) {
            zeldovichForLevel(0);
        } else {
            zeldovichForLevel(0);
            zeldovichForLevel(1);

            cerr << "Interpolating low-frequency information into zoom region...";
            cerr.flush();
            interpolateIntoLevel(1);
            cerr << "done." << endl;


            cerr << "Re-introducing high-k modes into low-res region...";
            recombineLevel0();
            zeldovichForLevel(0);
            cerr << "done." << endl;

        }
    }

    ///////////////////////////////////////////////
    // ZOOM particle management
    ///////////////////////////////////////////////

    void setInputMapper(std::string fname) {
      DummyIC<MyFloat> pseudoICs(this);
      auto dispatch = interpreter.specify_instance(pseudoICs);
      ifstream inf;
      inf.open(fname);


      if(!inf.is_open())
        throw std::runtime_error("Cannot open IC paramfile for relative_to command");
      cerr << "******** Running commands in" << fname << " to work out relationship ***********" << endl;
      dispatch.run_loop(inf);
      cerr << *(pseudoICs.pMapper) << endl;
      cerr << "******** Finished with" << fname << " ***********" << endl;
      pInputMapper = pseudoICs.pMapper;

    }

    std::shared_ptr<Grid<MyFloat>> getGridWithOutputOffset(int level=0) {
        auto gridForOutput = pGrid[level];

        if(xOffOutput!=0 || yOffOutput!=0 || zOffOutput!=0) {
            gridForOutput = std::make_shared<OffsetGrid<MyFloat>>(pGrid[0],xOffOutput,yOffOutput,zOffOutput);
        }
        return gridForOutput;
    }

    void initMapper() {

      if(pGrid.size()==0)
        return;

      // make a basic mapper for the base level grid
      pMapper = std::shared_ptr<ParticleMapper<MyFloat>>(new OneLevelParticleMapper<MyFloat>(
            getGridWithOutputOffset(0)
      ));


      if(pGrid.size()>=3) {
        // possible future enhancement, but for now...
        throw runtime_error("Don't know how to set up a mapper for more than one level of refinement");
      }

      if(pGrid.size()==2) {
        // it's a zoom!
        auto pMapperLevel1 = std::shared_ptr<ParticleMapper<MyFloat>>(new OneLevelParticleMapper<MyFloat>(getGridWithOutputOffset(1)));

        pMapper = std::shared_ptr<ParticleMapper<MyFloat>>(
            new TwoLevelParticleMapper<MyFloat>(pMapper, pMapperLevel1, zoomParticleArray, zoomfac*zoomfac*zoomfac));
      }

      if(Ob0>0) {

          // Add gas only to the deepest level. Pass the whole pGrid
          // vector if you want to add gas to every level.
          auto gasMapper = pMapper->addGas(Ob0/Om0,
                                          {pGrid.back()});


          // graft the gas particles onto the start of the map
          pMapper = std::make_shared<AddGasMapper<MyFloat>>(
              gasMapper.first, gasMapper.second, true);


      }

      // potentially resample the lowest-level DM grid. Again, this is theoretically
      // more flexible if you pass in other grid pointers.
      if(supersample>1)
          pMapper = pMapper->superOrSubSampleDM(supersample, {pGrid.back()},true);

      if(subsample>1)
        pMapper = pMapper->superOrSubSampleDM(subsample, {pGrid[0]}, false);

    }

    void clearAndDistributeParticleList() {

        if(pInputMapper!=nullptr) {
          pInputMapper->clearParticleList();
          pInputMapper->distributeParticleList(genericParticleArray);
        }
        else
        {
          pMapper->clearParticleList();
          pMapper->distributeParticleList(genericParticleArray);
        }
    }



    virtual void write() {

        zeldovich();

        cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
        cerr << (*pMapper);

        /*
        if (gadgetformat==2){
            SaveGadget2( (base+ ".gadget2").c_str(), nPartTotal,
            CreateGadget2Header<MyFloat>(nPartTotal, 0, a, zin, boxlen[0], Om0, Ol0, hubble),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3, Mass);
        }
        */
        if (gadgetformat==3){
            SaveGadget3( (base+ ".gadget3").c_str(), boxlen[0], Om0, Ol0, hubble, a, pMapper );
        }

        else if (gadgetformat==4)
            SaveTipsy(base+".tipsy", boxlen[0], Om0, Ol0, hubble, a, pMapper);

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

    int deepestLevelWithParticlesSelected() {
        if(pGrid.size()>1 && pGrid[1]->estimateParticleListSize()>0) return 1; else return 0;
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

        int level = deepestLevelWithParticlesSelected();

        MyFloat xa,ya,za, xb, yb, zb;

        std::vector<size_t> particleArray;
        pGrid[level]->gatherParticleList(particleArray);

        pGrid[level]->getCentroidLocation(particleArray[0],xa,ya,za);

        for(size_t i=0;i<particleArray.size();i++) {
            pGrid[level]->getCentroidLocation(particleArray[i],xb,yb,zb);
            x0+=get_wrapped_delta(xa,xb);
            y0+=get_wrapped_delta(ya,yb);
            z0+=get_wrapped_delta(za,zb);
        }
        x0/=particleArray.size();
        y0/=particleArray.size();
        z0/=particleArray.size();
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

        clearAndDistributeParticleList();
    }


    void cen_deriv4_alpha(long index, int direc, complex<MyFloat> *alpha, int level)

    {//4th order central difference

      MyFloat xp, yp, zp;//, zm2, zm1, zp2, zp1;
      pGrid[level]->getCentroidLocation(index, xp, yp, zp);
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

          ind_m1=pGrid[level]->findNextIndNoWrap(index, neg_step1);
          ind_p1=pGrid[level]->findNextIndNoWrap(index, step1);
          ind_m2=pGrid[level]->findNextIndNoWrap(ind_m1, neg_step1);
          ind_p2=pGrid[level]->findNextIndNoWrap(ind_p1, step1);

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

        std::vector<size_t> particleArray;
        pGrid[level]->gatherParticleList(particleArray);

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
            MyFloat w = 1.0/particleArray.size();
            for(size_t i=0;i<particleArray.size();i++) {
                rval[particleArray[i]]+=w;
            }

            fft(rval_k, rval, this->n[level], 1);
        }
        else if(strcasecmp(name,"phi")==0) {
            MyFloat w = 1.0/particleArray.size();
            for(size_t i=0;i<particleArray.size();i++) {
                rval[particleArray[i]]+=w;
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

            for(size_t i=0;i<particleArray.size();i++) {
                cen_deriv4_alpha(particleArray[i], direction, rval, level);
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

    void dumpID(string fname) {
        std::vector<size_t> results;
        pMapper->gatherParticleList(results);
        dumpBuffer(results, fname);
    }

    void centreParticle(long id) {
        pGrid[0]->getCentroidLocation(id,x0,y0,z0);
    }

    void selectNearest() {
        throw std::runtime_error("selectNearest not implemented");

    }

    void selectSphere(float radius) {
        MyFloat r2 = radius*radius;
        MyFloat delta_x, delta_y, delta_z, r2_i;
        MyFloat xp,yp,zp;

        genericParticleArray.clear();

        for_each_level(level) {
            std::vector<size_t> particleArray;

            for(size_t i=0;i<this->nPartLevel[level];i++) {
                pGrid[level]->getCentroidLocation(i,xp,yp,zp);
                delta_x = get_wrapped_delta(xp,x0);
                delta_y = get_wrapped_delta(yp,y0);
                delta_z = get_wrapped_delta(zp,z0);
                r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
                if(r2_i<r2)
                  particleArray.push_back(i);
            }

            pGrid[level]->clearParticleList();
            pGrid[level]->distributeParticleList(particleArray);

        }

    }


    void centreDenmax() {

        float den_max=-1000;
        long index_max=0;

        auto & pField_x = pGrid[0]->getFieldReal();
        for(size_t i=0;i<genericParticleArray.size();i++) {
            if(std::real(pField_x[genericParticleArray[i]])>den_max) {
                index_max=genericParticleArray[i];
                den_max = std::real(pField_x[genericParticleArray[i]]);
                pGrid[0]->getCentroidLocation(i,x0,y0,z0);
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
            pGrid[0]->getCentroidLocation(i,delta_x, delta_y, delta_z);
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
            val+=std::real(dot(vec, pGrid[level]->getFieldFourier().data(), this->nPartLevel[level]));
            cerr << val << " ";
            free(vec);
        }
        cerr << endl;
        cout << name << ": calculated value = " <<  val << endl;
    }

    virtual void constrain(string name, string type, float value) {
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


    virtual void fixConstraints() {
        if(pConstrainer==nullptr) {
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
                this->pGrid[level]->getFieldFourier()[i]=y1[j];
                j++;
            }
        }

        free(y1);

        cout << "Expected Delta chi^2=" << pConstrainer->get_delta_chi2() << endl;

        delete pConstrainer;
        pConstrainer=nullptr;
    }

    virtual void done() {
        if(pConstrainer!=nullptr) fixConstraints();
        write();
    }

    void reverse() {
        if(pConstrainer!=nullptr) {
            throw runtime_error("Can only reverse field direction after a 'fixconstraints' or 'done' command finailises the constraints");
        }

        for_each_level(level) {
            auto & field = pGrid[level]->getField();
            for(size_t i=0; i<this->nPartLevel[level]; i++)
                field[i]=-field[i];
        }
    }

    void reverseSmallK(MyFloat kmin) {
        if(pConstrainer!=nullptr) {
            throw runtime_error("Can only reverse field direction after a 'fixconstraints' or 'done' command finailises the constraints");
        }

        MyFloat k2min = kmin*kmin;


        for_each_level(level) {
	  MyFloat k2_g_min=std::numeric_limits<MyFloat>::max();
	  MyFloat k2_g_max=0.0;
	  size_t modes_reversed=0;
	  size_t tot_modes =pGrid[level]->size3;
	  auto & field = pGrid[level]->getFieldFourier();
	  const auto & grid = *(this->pGrid[level]);
	  MyFloat k2;
	  for(size_t i=0; i<this->nPartLevel[level]; i++) {
	    k2 = grid.getKSquared(i);
	    if(k2<k2min && k2!=0) {
	      field[i]=-field[i];
	      modes_reversed++;
	    }
	    if(k2<k2_g_min && k2!=0)
	      k2_g_min = k2;
	    if(k2>k2_g_max)
	      k2_g_max = k2;
	  }
	  cerr << "reverseSmallK: k reversal at " << sqrt(k2min) << "; grid was in range " << sqrt(k2_g_min) << " to " << sqrt (k2_g_max) << endl;
	  cerr << "               modes reversed = " << modes_reversed << " of " << tot_modes << endl;
        }
		
    }


};
#endif
