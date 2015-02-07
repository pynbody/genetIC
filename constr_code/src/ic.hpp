#include <string>
#include <tuple>
#include <gsl/gsl_rng.h> //link -lgsl and -lgslcblas at the very end
#include <gsl/gsl_randist.h> //for the gaussian (and other) distributions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;


template<typename MyFloat>
class IC {
protected:

    MyFloat Om0, Ol0, zin, a, sigma8, boxlen, dx;
    int out, n, gadgetformat, seed;

    int quoppa; // eh?
    double *kcamb, *Tcamb;
    long nPartTotal;
    string incamb, indir, inname, base;

    bool whiteNoiseFourier;

    complex<MyFloat> *pField_k;
    complex<MyFloat> *P;

    Grid<MyFloat> *pGrid;
    int *part_arr;
    int n_part_arr;
    MyFloat x0, y0, z0;
    complex<MyFloat> *pField_x;
    MultiConstrainedField<MyFloat> *pConstrainer;



public:
    IC() {
        pGrid=NULL;
        pField_x=NULL;
        pConstrainer=NULL;
        part_arr=NULL;
        n_part_arr=0;
        whiteNoiseFourier=false;
    }

    ~IC() {
        if(pGrid!=NULL) delete pGrid;

    }

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
        boxlen = in;
        dx = boxlen/n;
    }

    void setZ0(MyFloat in) {
        zin = in;
        a=1./(zin+1.);
    }

    void setn(int in) {
        n = in;
        dx = boxlen/n;
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
        prepare();
    }

    string make_base( string basename, int n, MyFloat box, MyFloat zin){
        ostringstream nult;
        if(inname.size()==0) {
            nult << basename<<"IC_iter_" << floatinfo<MyFloat>::name << "_z"<<zin<<"_"<<n<<"_L" << box;

        } else {
            nult << basename << "/" << inname;
        }
        return nult.str();
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

    long kgridToIndex(int k1, int k2, int k3) const {
        long ii,jj,ll,idk;
        if(k1<0) ii=k1+n; else ii=k1;
        if(k2<0) jj=k2+n; else jj=k2;
        if(k3<0) ll=k3+n; else ll=k3;
        idk=(ii*(n)+jj)*(n)+ll;
        return idk;
    }

    void drawOneFourierMode(gsl_rng *r, int k1, int k2, int k3, MyFloat norm) {
        long id_k, id_negk;
        id_k = kgridToIndex(k1,k2,k3);
        id_negk = kgridToIndex(-k1,-k2,-k3);
        pField_k[id_k]=std::complex<MyFloat>(norm*gsl_ran_gaussian_ziggurat(r,1.),norm*gsl_ran_gaussian_ziggurat(r,1.));

        // reality condition:
        pField_k[id_negk]=std::conj(pField_k[id_k]);
    }

    void drawRandom() {
        gsl_rng * r;
        const gsl_rng_type * T; //generator type variable

        T = gsl_rng_ranlxs2; // shouldn't this be gsl_rng_ranlxd2 for MyFloat = double?
        r = gsl_rng_alloc (T); //this allocates memory for the generator with type T
        gsl_rng_set(r,seed);

        pField_k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        MyFloat sigma=sqrt((MyFloat)(nPartTotal));

        if(!whiteNoiseFourier) {
            // original implementation


            complex<MyFloat> *rnd=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

            cerr << "Drawing random numbers..."<< endl;

            long i;

            for(i=0;i<nPartTotal;i++){rnd[i]=gsl_ran_gaussian_ziggurat(r,1.)*sigma;}// cout<< "rnd "<< rnd[i] << endl;}

            gsl_rng_free (r);

            cout<< "First FFT..." <<endl;
            fft_r(pField_k, rnd, n, 1);

            free(rnd);

            pField_k[0]=complex<MyFloat>(0.,0.); //assure mean==0

        } else {

            // do it in fourier space, in order of increasing |k|, so that
            // resolution can be scaled by factors of 2 and we still get
            // the 'same' field

            cerr << "Drawing random numbers in fourier space..."<< endl;
            MyFloat kk=0.;
            int ks,k1,k2,k3;

            sigma/=sqrt(2.0);

            // Do it in square k-shells
            for(ks=0; ks<n/2;ks++)
            {
                for(k1=-ks; k1<ks; k1++) {
                    for(k2=-ks; k2<ks; k2++) {
                        drawOneFourierMode(r,ks,k1,k2,sigma);
                        drawOneFourierMode(r,k1,ks,k2,sigma);
                        drawOneFourierMode(r,k1,k2,ks,sigma);

                    }
                }
            }

        }


    }

    void applyPowerSpec() {
        //scale white-noise delta with initial PS


        MyFloat grwfac;
        MyFloat ns=0.96; //TODO make this an input parameter (long term: run CAMB with input parameters to ensure consistency?)

        //growth factor normalised to 1 today:
        grwfac=D(a, Om0, Ol0)/D(1., Om0, Ol0);

        cout<< "Growth factor " << grwfac << endl;
        MyFloat sg8;

        sg8=sig(8., kcamb, Tcamb, ns, boxlen, n, quoppa);
        std::cout <<"Sigma_8 "<< sg8 << std::endl;

        MyFloat kw = 2.*M_PI/(MyFloat)boxlen;
        MyFloat amp=(sigma8/sg8)*(sigma8/sg8)*grwfac*grwfac; //norm. for sigma8 and linear growth factor
        MyFloat norm=kw*kw*kw/powf(2.*M_PI,3.); //since kw=2pi/L, this is just 1/V_box

        P=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));


        std::cout<<"Interpolation: kmin: "<< kw <<" Mpc/h, kmax: "<< (MyFloat)( kw*n/2.*sqrt(3.)) <<" Mpc/h"<<std::endl;

        MyFloat norm_amp=norm*amp;

        brute_interpol_new(n, kcamb, Tcamb,  quoppa, kw, ns, norm_amp, pField_k, pField_k, P);

        //assert(abs(real(norm_iter)) >1e-12); //norm!=0
        //assert(abs(imag(norm_iter)) <1e-12); //because there shouldn't be an imaginary part since we assume d is purely real
        //TODO think about this (it fails more often than it should?)

        cout<<"Transfer applied!"<<endl;
        cout<< "Power spectrum sample: " << P[0] << " " << P[1] <<" " << P[nPartTotal-1] <<endl;

        std::cout<<"Initial chi^2 (white noise, fourier space) = " << chi2(pField_k,P,nPartTotal) << std::endl;
        powsp_noJing(n, pField_k, (base+ ".ps").c_str(), boxlen);
    }

    std::tuple<MyFloat*, MyFloat*, MyFloat*, MyFloat*, MyFloat*, MyFloat*> zeldovich() {
        //Zeldovich approx.

        complex<MyFloat>* psift1k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        complex<MyFloat>* psift2k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        complex<MyFloat>* psift3k=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

        int iix, iiy, iiz;
        long idx;
        MyFloat kfft;
        MyFloat kw = 2.*M_PI/(MyFloat)boxlen;

        for(int ix=0; ix<n;ix++){
            for(int iy=0;iy<n;iy++){
                for(int iz=0;iz<n;iz++){


                    idx = (ix*n+iy)*(n)+iz;

                    if( ix>n/2 ) iix = ix - n; else iix = ix;
                    if( iy>n/2 ) iiy = iy - n; else iiy = iy;
                    if( iz>n/2 ) iiz = iz - n; else iiz = iz;

                    kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);

                    psift1k[idx].real(-pField_k[idx].imag()/(MyFloat)(kfft*kfft)*iix/kw);
                    psift1k[idx].imag(pField_k[idx].real()/(MyFloat)(kfft*kfft)*iix/kw);
                    psift2k[idx].real(-pField_k[idx].imag()/(MyFloat)(kfft*kfft)*iiy/kw);
                    psift2k[idx].imag(pField_k[idx].real()/(MyFloat)(kfft*kfft)*iiy/kw);
                    psift3k[idx].real(-pField_k[idx].imag()/(MyFloat)(kfft*kfft)*iiz/kw);
                    psift3k[idx].imag(pField_k[idx].real()/(MyFloat)(kfft*kfft)*iiz/kw);
                }
            }
        }


        /* //Optional stuff:

        complex<MyFloat> *delta_real=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        delta_real=fft_r(delta_real,pField_k,n,-1);

        MyFloat *potr=(MyFloat*)calloc(n*n*n,sizeof(MyFloat));
        for(i=0;i<n*n*n;i++) {
            potr[i]=(MyFloat)(pot[i].real());
            pot[i].imag(0.);
        }

        //optional: save potential
        //string phases_out=(base+ "_phases.hdf5").c_str();
        //save_phases(potk, potr, delta_real, n, phases_out.c_str());

        MyFloat exp=0.;
        for(i=0;i<n*n*n;i++) exp+=(MyFloat)(pot[i].real()*pot[i].real());

        exp/=(MyFloat)(n*n*n);

        free(pot);
        free(potk);
        free(potr);

        free(delta_real);

        free(pField_k);
        */

        psift1k[0]=complex<MyFloat>(0.,0.);
        psift2k[0]=complex<MyFloat>(0.,0.);
        psift3k[0]=complex<MyFloat>(0.,0.);

        complex<MyFloat>* psift1=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        complex<MyFloat>* psift2=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));
        complex<MyFloat>* psift3=(complex<MyFloat>*)calloc(nPartTotal,sizeof(complex<MyFloat>));

        psift1=fft_r(psift1,psift1k,n,-1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
        psift2=fft_r(psift2,psift2k,n,-1); //same
        psift3=fft_r(psift3,psift3k,n,-1); //same

        free(psift1k);
        free(psift2k);
        free(psift3k);

        MyFloat gr=boxlen/(MyFloat)n;
        cout<< "Grid cell size: "<< gr <<" Mpc/h"<<endl;

        MyFloat *Vel1=(MyFloat*)calloc(nPartTotal,sizeof(MyFloat));
        MyFloat *Vel2=(MyFloat*)calloc(nPartTotal,sizeof(MyFloat));
        MyFloat *Vel3=(MyFloat*)calloc(nPartTotal,sizeof(MyFloat));
        MyFloat *Pos1=(MyFloat*)calloc(nPartTotal,sizeof(MyFloat));
        MyFloat *Pos2=(MyFloat*)calloc(nPartTotal,sizeof(MyFloat));
        MyFloat *Pos3=(MyFloat*)calloc(nPartTotal,sizeof(MyFloat));

        MyFloat hfac=1.*100.*sqrt(Om0/a/a/a+Ol0)*sqrt(a);
        //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
        //TODO: hardcoded value of f=1 is inaccurate, but fomega currently gives wrong nults


        MyFloat mean1=0.,mean2=0.,mean3=0.;
        MyFloat Mean1=0., Mean2=0., Mean3=0.;
        cout<< "Applying ZA & PBC... "<<endl;
        //apply ZA:
        for(int ix=0;ix<n;ix++) {
            for(int iy=0;iy<n;iy++) {
                for(int iz=0;iz<n;iz++) {

                    idx = (ix*n+iy)*(n)+iz;

                    Vel1[idx] = psift1[idx].real()*hfac; //physical units
                    Vel2[idx] = psift2[idx].real()*hfac;
                    Vel3[idx] = psift3[idx].real()*hfac;

                    //position in "grid coordinates": Pos e [0,N-1]
                    Pos1[idx] = psift1[idx].real()*n/boxlen+ix;
                    Pos2[idx] = psift2[idx].real()*n/boxlen+iy;
                    Pos3[idx] = psift3[idx].real()*n/boxlen+iz;

                    mean1+=abs(psift1[idx].real()*n/boxlen);
                    mean2+=abs(psift2[idx].real()*n/boxlen);
                    mean3+=abs(psift3[idx].real()*n/boxlen);

                    //enforce periodic boundary conditions
                    if(Pos1[idx]<0.) Pos1[idx]+=(MyFloat)n;
                    else if(Pos1[idx]>(MyFloat)n) Pos1[idx]-=(MyFloat)n;
                    if(Pos2[idx]<0.) Pos2[idx]+=(MyFloat)n;
                    else if(Pos2[idx]>(MyFloat)n) Pos2[idx]-=(MyFloat)n;
                    if(Pos3[idx]<0.) Pos3[idx]+=(MyFloat)n;
                    else if(Pos3[idx]>(MyFloat)n) Pos3[idx]-=(MyFloat)n;

                    //ncale to physical coordinates
                    Pos1[idx] *= (boxlen/(MyFloat)n);
                    Pos2[idx] *= (boxlen/(MyFloat)n);
                    Pos3[idx] *= (boxlen/(MyFloat)n);

                    Mean1+=Pos1[idx];
                    Mean2+=Pos2[idx];
                    Mean3+=Pos3[idx];

                }
            }
        }


        cout<< "Box/2="<< boxlen/2.<< " Mpc/h, Mean position x,y,z: "<< Mean1/(MyFloat(n*n*n))<<" "<< Mean2/(MyFloat(n*n*n))<<" "<<Mean3/(MyFloat(n*n*n))<< " Mpc/h"<<  endl;

        free(psift1);
        free(psift2);
        free(psift3);

        return make_tuple(Pos1,Pos2,Pos3,Vel1,Vel2,Vel3);


    }

    void write() {

        MyFloat *Pos1, *Pos2, *Pos3, *Vel1, *Vel2, *Vel3;

        tie(Pos1,Pos2,Pos3,Vel1,Vel2,Vel3) = zeldovich();

        MyFloat pmass=27.78*Om0*powf(boxlen/(MyFloat)(n),3.0); // in 10^10 M_sol
        cout<< "Particle mass: " <<pmass <<" [10^10 M_sun]"<<endl;

        if (gadgetformat==2){
            SaveGadget2( (base+ "_gadget2.dat").c_str(), n,
            CreateGadget2Header<MyFloat>(nPartTotal, pmass, a, zin, boxlen, Om0, Ol0),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
        }

        if (gadgetformat==3){
            SaveGadget3( (base+ "_gadget3.dat").c_str(), n,
            CreateGadget3Header<MyFloat>(nPartTotal, pmass, a, zin, boxlen, Om0, Ol0),
            Pos1, Vel1, Pos2, Vel2, Pos3, Vel3);
        }

        free(Pos1);
        free(Vel1);
        free(Pos2);
        free(Vel2);
        free(Pos3);
        free(Vel3);
    }


    virtual void prepare() {
        nPartTotal = ((long)n*n)*n;

        base = make_base(indir,n,boxlen,zin);

        readCamb();
        drawRandom();
        applyPowerSpec();

        pGrid = new Grid<MyFloat>(n);

        UnderlyingField<MyFloat> *pField =
            new UnderlyingField<MyFloat>(this->P, this->pField_k, this->nPartTotal);

        pConstrainer = new MultiConstrainedField<MyFloat>(pField, this->nPartTotal);
    }

    ///////////////////////
    // CONSTRAING CODE //
    ///////////////////////



protected:


    MyFloat get_wrapped_delta(MyFloat x0, MyFloat x1) {
        MyFloat result = x0-x1;
        if(result>this->boxlen/2) {
            result-=this->boxlen;
        }
        if(result<-this->boxlen/2) {
            result+=this->boxlen;
        }
        return result;
    }


    void getCentre() {

        x0 = 0; y0 = 0; z0 =0;

        MyFloat xa=this->pGrid->cells[part_arr[0]].coords[0]*this->dx;
        MyFloat ya=this->pGrid->cells[part_arr[0]].coords[1]*this->dx;
        MyFloat za=this->pGrid->cells[part_arr[0]].coords[2]*this->dx;

        for(long i=0;i<n_part_arr;i++) {
            x0+=get_wrapped_delta(this->pGrid->cells[part_arr[i]].coords[0]*this->dx,xa);
            y0+=get_wrapped_delta(this->pGrid->cells[part_arr[i]].coords[1]*this->dx,ya);
            z0+=get_wrapped_delta(this->pGrid->cells[part_arr[i]].coords[2]*this->dx,za);
        }
        x0/=n_part_arr;
        y0/=n_part_arr;
        z0/=n_part_arr;
        x0+=xa;
        y0+=ya;
        z0+=za;
    }




    void AllocAndGetBuffer_int(const char* IDfile, bool append=false) {
        // count lines, allocate space, then read

        FILE *f;
        int i=0;
        int r=0;
        int c;

        f = fopen(IDfile, "r");
        if (f == NULL) throw std::runtime_error("File not found");

        while ( (c=fgetc(f)) != EOF ) {
            if ( c == '\n' )
                    i++;
        }

        fclose(f);

        cerr << "File " <<IDfile << " has " << i << " lines" << endl;

        if(append && part_arr==NULL) {
            cerr << "Can't append - no particles were loaded" << endl;
            append=false;
        }
        int final_size = i;

        if(append) final_size+=n_part_arr;
        int *arr = (int*)calloc(final_size,sizeof(int));

        GetBuffer_int(arr, IDfile, i);

        // copy old particles into new buffer
        if (append>0)
          cerr << "Also keeping " << n_part_arr << " existing particles" << endl;
        for(int j=0; j<n_part_arr; j++)
          arr[i+j] = part_arr[j];

        if(part_arr!=NULL)
          free(part_arr);

        part_arr = arr;
        n_part_arr = final_size;

        getCentre();

    }


    void cen_deriv4_alpha(long index, int direc, complex<MyFloat> *alpha)

    {//4th order central difference

      MyFloat x0, y0, z0;//, zm2, zm1, zp2, zp1;
      x0=get_wrapped_delta(dx*pGrid->cells[index].coords[0],x0);
      y0=get_wrapped_delta(dx*pGrid->cells[index].coords[1],y0);
      z0=get_wrapped_delta(dx*pGrid->cells[index].coords[2],z0);

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

          ind_m1=pGrid->find_next_ind(index, neg_step1);
          ind_p1=pGrid->find_next_ind(index, step1);


          ind_m2=pGrid->find_next_ind(ind_m1, neg_step1);
          ind_p2=pGrid->find_next_ind(ind_p1, step1);

          MyFloat a=-1./12./dx, b=2./3./dx;  //the signs here so that L ~ - Nabla Phi

          alpha[ind_m2]+=(c[di]*a);
          alpha[ind_m1]+=(c[di]*b);
          alpha[ind_p1]+=(-c[di]*b);
          alpha[ind_p2]+=(-c[di]*a);
       }

    }

    complex<MyFloat> *calcConstraintVector(string name_in) {
        const char* name = name_in.c_str();

        complex<MyFloat> *rval=(complex<MyFloat>*)calloc(this->nPartTotal,sizeof(complex<MyFloat>));
        complex<MyFloat> *rval_k=(complex<MyFloat>*)calloc(this->nPartTotal,sizeof(complex<MyFloat>));

        if(strcasecmp(name,"overdensity")==0) {
            MyFloat w = 1.0/n_part_arr;
            for(long i=0;i<n_part_arr;i++) {
                rval[part_arr[i]]=w;
            }

            fft_r(rval_k, rval, this->n, 1);
        }
        else if(strcasecmp(name,"phi")==0) {
            MyFloat w = 1.0/n_part_arr;
            for(long i=0;i<n_part_arr;i++) {
                rval[part_arr[i]]=w;
            }
            complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(this->nPartTotal,sizeof(complex<MyFloat>));
            fft_r(rval_kX, rval, this->n, 1);
            poiss(rval_k, rval_kX, this->n, this->boxlen, this->a, this->Om0);
            free(rval_kX);
        }
        else if(name[0]=='L' || name[0]=='l') {
            // angular momentum
            int direction=atoi(&name[1]);

            cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;

            for(long i=0;i<n_part_arr;i++) {
                cen_deriv4_alpha(part_arr[i], direction, rval);
            }
            complex<MyFloat> *rval_kX=(complex<MyFloat>*)calloc(this->nPartTotal,sizeof(complex<MyFloat>));
            fft_r(rval_kX, rval, this->n, 1);
            // The constraint as derived is on the potential. By considering
            // unitarity of FT, we can FT the constraint to get the constraint
            // on the density.
            poiss(rval_k, rval_kX, this->n, this->boxlen, this->a, this->Om0);
            free(rval_kX);
        } else {
            cout << "  -> UNKNOWN constraint vector type, returning zeros [this is bad]" << endl;

        }

        if(pField_x!=NULL) {
            cout << " dot in real space = " << std::real(dot(pField_x, rval, this->nPartTotal)) << endl;
        }

        free(rval);
        return rval_k;

    }

    void ensureRealDelta() {
        if(pField_x==NULL) {
            pField_x = (complex<MyFloat>*)calloc(this->nPartTotal,sizeof(complex<MyFloat>));
            fft_r(pField_x,this->pField_k,this->n,-1);
        }
    }



public:


    void loadID(string fname) {
        AllocAndGetBuffer_int(fname.c_str());
    }

    void appendID(string fname) {
        AllocAndGetBuffer_int(fname.c_str(), true);
    }

    void centreParticle(long id) {
        x0 = this->pGrid->cells[id].coords[0]*this->dx;
        y0 = this->pGrid->cells[id].coords[1]*this->dx;
        z0 = this->pGrid->cells[id].coords[2]*this->dx;
    }

    void selectSphere(float radius) {
        float r2 = radius*radius;
        float delta_x, delta_y, delta_z, r2_i;
        int n=0;
        for(long i=0;i<this->nPartTotal;i++) {
            delta_x = get_wrapped_delta(this->pGrid->cells[i].coords[0]*this->dx,x0);
            delta_y = get_wrapped_delta(this->pGrid->cells[i].coords[1]*this->dx,y0);
            delta_z = get_wrapped_delta(this->pGrid->cells[i].coords[2]*this->dx,z0);
            r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
            if(r2_i<r2)
                n++;
        }
        cerr << "Selecting " << n << " particles..." << endl;
        if(part_arr!=NULL)
            free(part_arr);

        part_arr = (int*)calloc(n,sizeof(int));
        n=0;
        for(long i=0;i<this->nPartTotal;i++) {
            delta_x = get_wrapped_delta(this->pGrid->cells[i].coords[0]*this->dx,x0);
            delta_y = get_wrapped_delta(this->pGrid->cells[i].coords[1]*this->dx,y0);
            delta_z = get_wrapped_delta(this->pGrid->cells[i].coords[2]*this->dx,z0);
            r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
            if(r2_i<r2) part_arr[n++]=i;
        }
        n_part_arr = n;
    }


    void centreDenmax() {

        float den_max=-1000;
        long index_max=0;
        ensureRealDelta();

        for(long i=0;i<n_part_arr;i++) {
            if(std::real(pField_x[part_arr[i]])>den_max) {
                index_max=part_arr[i];
                den_max = std::real(pField_x[part_arr[i]]);
                x0 = this->pGrid->cells[part_arr[i]].coords[0]*this->dx;
                y0 = this->pGrid->cells[part_arr[i]].coords[1]*this->dx;
                z0 = this->pGrid->cells[part_arr[i]].coords[2]*this->dx;
            }
        }

        cerr << "Denmax = " << den_max <<", index=" << index_max << " coords=" << x0 << " " << y0 << " " << z0 << endl;

    }

    void reorderBuffer() {

        cout << "Reordering buffer radially..." << endl;
        cout << " [taking centre = " << x0 << " " << y0 << " " <<z0 << "]" << endl;

        MyFloat r2[n_part_arr];
        MyFloat delta_x, delta_y, delta_z;
        std::vector<size_t> index(n_part_arr);

        for(int i=0;i<n_part_arr;i++) {
            delta_x = get_wrapped_delta(this->pGrid->cells[part_arr[i]].coords[0]*this->dx,x0);
            delta_y = get_wrapped_delta(this->pGrid->cells[part_arr[i]].coords[1]*this->dx,y0);
            delta_z = get_wrapped_delta(this->pGrid->cells[part_arr[i]].coords[2]*this->dx,z0);
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

    void truncateBuffer(float x) {
        if(x<0 ) {
            cerr << "Truncate command takes a fraction between 0 and 1 or a number>=1" << endl;
            exit(0);
        }
        if(x<1)
            n_part_arr = ((int)(n_part_arr*x));
        else
            n_part_arr = ((int)x);
    }

    void calculate(string name) {
        complex<MyFloat> *vec = calcConstraintVector(name);
        cout << name << ": calculated value = " << dot(vec, this->pField_k, this->nPartTotal) << endl;
        free(vec);
    }

    void constrain(string name, string type, float value) {
        if(pConstrainer==NULL) {
            throw runtime_error("No constraint information is available. Is your done command too early, or repeated?");
        }

        bool relative;
        if (strcasecmp(type.c_str(),"relative")==0) {
            relative=true;
        } else if (strcasecmp(type.c_str(),"absolute")!=0) {
            cerr << "Constraints must state either relative or absolute" << endl;
            exit(0);
        }

        std::complex<MyFloat> constraint = value;
        complex<MyFloat> *vec = calcConstraintVector(name);
        std::complex<MyFloat> initv = dot(vec, this->pField_k, this->nPartTotal);

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
        complex<MyFloat> *y1=(complex<MyFloat>*)calloc(this->nPartTotal,sizeof(complex<MyFloat>));
        pConstrainer->prepare();
        pConstrainer->get_realization(y1);

        for(long i=0; i<this->nPartTotal; i++)
            this->pField_k[i]=y1[i];

        free(y1);

        cout << "Expected Delta chi^2=" << pConstrainer->get_delta_chi2() << endl;

        delete pConstrainer;
        pConstrainer=NULL;

    }


};
