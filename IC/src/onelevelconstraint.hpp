//
// Created by Andrew Pontzen on 15/06/15.
//

#ifndef IC_CONSTRAINT_HPP
#define IC_CONSTRAINT_HPP

#include <string>
#include <complex>
#include <vector>
#include <map>
#include "grid.hpp"
#include "cosmo.hpp"

template<typename MyFloat>
class ConstraintCalculator {
protected:
    const Grid<MyFloat> &grid;
    std::vector<complex<MyFloat>> &output;
    std::vector<size_t> particleArray;
    const CosmologicalParameters<MyFloat> &cosmology;
    MyFloat x0=0., y0=0., z0=0.; //centre coordinates needed for angmom; set in getCentre()

    void cen_deriv4_alpha(long index, int direc, MyFloat x0, MyFloat y0, MyFloat z0)

    {   //4th order central difference

        MyFloat xp=0., yp=0., zp=0.;
        //grid.getCentroidLocation(index, xp, yp, zp);


        int xi, yi, zi;
        grid.getCoordinates(index,xi,yi,zi);

        xp=MyFloat(xi)*grid.dx;
        yp=MyFloat(yi)*grid.dx;
        zp=MyFloat(zi)*grid.dx;
        //xp=0;
        //yp=0;
        //zp=0.195312;

        //cerr<< "index" << index<< endl;

        //cerr << "centroid " <<  xp << " "<< yp << " " << zp << endl;

        xp=grid.getWrappedDelta(xp,x0);
        yp=grid.getWrappedDelta(yp,y0);
        zp=grid.getWrappedDelta(zp,z0);

        //cerr << "centroid after wrap" <<  xp << " "<< yp << " " << zp << endl;

        MyFloat c[3]={0,0,0};
        if(direc==0){c[2]=yp; c[1]=-zp;}
        else if(direc==1){c[0]=zp; c[2]=-xp;}
        else if(direc==2){c[1]=xp; c[0]=-yp;}
        else if(direc==3){
            MyFloat rp = std::sqrt((xp*xp)+(yp*yp)+(zp*zp));
            if(rp!=0) {
                c[0]=xp/rp;
                c[1]=yp/rp;
                c[2]=zp/rp;
            }
        } // radial velocity

        //do the epsilon
        // MyFloat c[3]={0,0,0};
        // if(direc==0){c[2]=y0; c[1]=-z0;}
        // else if(direc==1){c[0]=z0; c[2]=-x0;}
        // else if(direc==2){c[1]=x0; c[0]=-y0;}
        // else if(direc==3){
        //     MyFloat r0 = std::sqrt((x0*x0)+(y0*y0)+(z0*z0));
        //     if(r0!=0) {
        //         c[0]=x0/r0;
        //         c[1]=y0/r0;
        //         c[2]=z0/r0;
        //     }
        // } // radial velocity

        else{cerr<< "Wrong value for parameter 'direc' in function 'cen_deriv4_alpha'."<< endl; exit(1);}


        for(int di=0; di<3; di++) {
            long ind_p1, ind_m1, ind_p2, ind_m2;
            //first step in rho direction
            int step1[3]={0,0,0};
            int neg_step1[3]={0,0,0};
            step1[di]=1;
            neg_step1[di]=-1;

            // N.B. can't wrap - might be on subgrid

            //ind_m1=grid.findNextIndNoWrap(index, neg_step1);
            //ind_p1=grid.findNextIndNoWrap(index, step1);
            //ind_m2=grid.findNextIndNoWrap(ind_m1, neg_step1);
            //ind_p2=grid.findNextIndNoWrap(ind_p1, step1);


            ind_m1=grid.findNextInd(index, neg_step1);
            ind_p1=grid.findNextInd(index, step1);
            ind_m2=grid.findNextInd(ind_m1, neg_step1);
            ind_p2=grid.findNextInd(ind_p1, step1);


            MyFloat a=-1./12./grid.dx, b=2./3./grid.dx;  //the signs here so that L ~ - Nabla Phi

            //cerr << "direc " << di << " " << ind_m1 << " " << ind_m2 << " " << ind_p1 << " " << ind_p2 << endl;
            //cerr << "c[direc] " << c[di]/50.*256 <<  " c[direc]*(a,b): (" << c[di]*a/50.*256 << ", " << c[di]*b/50.*256 << ")" << endl;
            //cerr << endl;

            output[ind_m2]+=(c[di]*a);
            output[ind_m1]+=(c[di]*b);
            output[ind_p1]+=(-c[di]*b);
            output[ind_p2]+=(-c[di]*a);
        }

    }

public:
    ConstraintCalculator(const Grid<MyFloat>& grid, std::vector<complex<MyFloat>> & output,
                         const CosmologicalParameters<MyFloat> &cosmology)
            : grid(grid), output(output), cosmology(cosmology)
    {
        grid.gatherParticleList(particleArray);
        output.resize(grid.size3);
        getCentre();    
    }    

    void getCentre() {
        
        MyFloat xa,ya,za, xb, yb, zb, x0=0., y0=0., z0=0.;

        std::vector<size_t> particleArray;
        grid.gatherParticleList(particleArray);

        grid.getCentroidLocation(particleArray[0],xa,ya,za);
        //xp=grid.getWrappedDelta(xp,x0);


        cout << "xa  "<< xa << " " << ya << " " << za << endl;

        for(size_t i=0;i<particleArray.size();i++) {
            grid.getCentroidLocation(particleArray[i],xb,yb,zb);
            x0+=grid.getWrappedDelta(xa,xb);
            y0+=grid.getWrappedDelta(ya,yb);
            z0+=grid.getWrappedDelta(za,zb);
        }

        x0/=particleArray.size();
        y0/=particleArray.size();
        z0/=particleArray.size();
        x0+=xa;
        y0+=ya;
        z0+=za;

        cerr << "Inside CC: Centre of region is " << x0 << " " << y0 << " " << z0 << endl;

        //x0=0;
        //y0=0;
        //z0=grid.dx;

        x0=44.3359375;
        y0=33.3984375;
        z0=35.3515625;

        cerr<< "setting to " << x0 << " " << y0 << " "<< z0 << endl;

        this->x0=x0;
        this->y0=y0;
        this->z0=z0;

    }

    void overdensity() {

        MyFloat w = 1.0 / particleArray.size();

        for (size_t i = 0; i < grid.size3; ++i) {
            output[i] = 0;
        }

        for (size_t i = 0; i < particleArray.size(); i++) {
            output[particleArray[i]] += w;
        }

        fft(output.data(), output.data(), grid.size, 1);
    }

    void phi() {

        MyFloat w = 1.0 / particleArray.size();


        for (size_t i = 0; i < grid.size3; ++i) {
            output[i] = 0;
        }

        for (size_t i = 0; i < particleArray.size(); i++) {
            output[particleArray[i]] += w;
        }

        fft(output.data(), output.data(), grid.size, 1);
        poiss(output.data(), output.data(), grid.size, grid.boxsize, cosmology.scalefactor, cosmology.OmegaM0);
    }

//old angular momentum class, temporary fix below:
    void angmom(int direction) {
        //throw std::runtime_error("Angmom cannot run - x0,y0,z0 coordinates not yet passed to new constraint class");
        //cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;

        for(size_t i=0;i<particleArray.size();i++) {
            cen_deriv4_alpha(particleArray[i], direction, x0, y0, z0);
        }

        fft(output.data(), output.data(), grid.size, 1);
        // The constraint as derived is on the potential. By considering
        // unitarity of FT, we can FT the constraint to get the constraint
        // on the density.
        poiss(output.data(), output.data(), grid.size, grid.boxsize,
              cosmology.scalefactor, cosmology.OmegaM0);

    }

    void angmom0() {

        MyFloat x0, y0, z0;
        x0=this->x0;
        y0=this->y0;
        z0=this->z0;

        
        cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;

        // x0=44.3359375;
        // y0=33.3984375;
        // z0=35.3515625;
        // cerr << "Angmom centre is SET TO " <<x0 << " " <<y0 << " " << z0 << endl;
        
        for(size_t i=0;i<particleArray.size();i++) {
            cen_deriv4_alpha(particleArray[i], 0, x0, y0, z0);
        }

         //output[5]=complex<MyFloat> (0.666667,0);
         //output[7]=complex<MyFloat> (-0.666667,0);
         //output[13]=complex<MyFloat> (-0.666667,0);
         //output[15]=complex<MyFloat> (0.666667,0);

        int ii;
        int sum=0;
        for(ii=0; ii< grid.size3; ii++){
            if ( (output[ii]!=complex<MyFloat> (0.,0.)) && (sum<10)  ){
                cerr<< " output angmom 0 before "<< ii<< " " << output[ii] <<endl; sum+=1;} 
        }

        fft(output.data(), output.data(), grid.size, 1);

        //int ii;
        sum=0;
        for(ii=0; ii< grid.size3; ii++){
            if ( (output[ii]!=complex<MyFloat> (0.,0.)) && (sum<10)  ){
                cerr<< " output angmom 0 after "<< ii<< " " << output[ii] <<endl; sum+=1;} 
        }
        // The constraint as derived is on the potential. By considering
        // unitarity of FT, we can FT the constraint to get the constraint
        // on the density.
        poiss(output.data(), output.data(), grid.size, grid.boxsize,
              cosmology.scalefactor, cosmology.OmegaM0);


        sum=0;
        for(ii=0; ii< grid.size3; ii++){
            if ( (output[ii]!=complex<MyFloat> (0.,0.)) && (sum<10)  ){
                cerr<< " output angmom 0 end "<< ii<< " " << output[ii] <<endl; sum+=1;}
        } 

    }

    void angmom1() {
        
        MyFloat x0, y0, z0;
        x0=this->x0;
        y0=this->y0;
        z0=this->z0;

        cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;

        // x0=44.3359375;
        // y0=33.3984375;
        // z0=35.3515625;
        // cerr << "Angmom centre is SET TO " <<x0 << " " <<y0 << " " << z0 << endl;

        //for(size_t i=0;i<10;i++) {cout<< "part arr, derivative1: " <<  particleArray[i] << " "<< endl;}//cen_deriv4_alpha(particleArray[i], 1)<< endl;}

        for(size_t i=0;i<particleArray.size();i++) {
            cen_deriv4_alpha(particleArray[i], 1, x0, y0, z0);
        }

        fft(output.data(), output.data(), grid.size, 1);
        // The constraint as derived is on the potential. By considering
        // unitarity of FT, we can FT the constraint to get the constraint
        // on the density.
        poiss(output.data(), output.data(), grid.size, grid.boxsize,
              cosmology.scalefactor, cosmology.OmegaM0);

    }

    void angmom2() {
        
        MyFloat x0, y0, z0;
        x0=this->x0;
        y0=this->y0;
        z0=this->z0;

        cerr << "Angmom centre is " <<x0 << " " <<y0 << " " << z0 << endl;


        // x0=44.3359375;
        // y0=33.3984375;
        // z0=35.3515625;
        // cerr << "Angmom centre is SET TO " <<x0 << " " <<y0 << " " << z0 << endl;

        for(size_t i=0;i<particleArray.size();i++) {
            cen_deriv4_alpha(particleArray[i], 2, x0, y0, z0);
        }

        fft(output.data(), output.data(), grid.size, 1);
        // The constraint as derived is on the potential. By considering
        // unitarity of FT, we can FT the constraint to get the constraint
        // on the density.
        poiss(output.data(), output.data(), grid.size, grid.boxsize,
              cosmology.scalefactor, cosmology.OmegaM0);

    }
};



template<typename MyFloat>
void calcConstraint(const std::string &name, const Grid<MyFloat>& grid,
                    const CosmologicalParameters<MyFloat> &cosmology,
                    std::vector<complex<MyFloat>> &output)
{
    typedef ConstraintCalculator<MyFloat> CC;

    CC calc(grid,output,cosmology);

    
    InstanceDispatch<CC,void> dispatch(calc);    

    dispatch.add_class_route("overdensity",&CC::overdensity);
    dispatch.add_class_route("phi",&CC::phi);
    //temporary and dirty fix to calculate angular momentum (by Nina):
    dispatch.add_class_route("L0",&CC::angmom0);
    dispatch.add_class_route("L1",&CC::angmom1);
    dispatch.add_class_route("L2",&CC::angmom2);

    dispatch.run(name);

    int ii;
    int sum=0;
    for(ii=0; ii< grid.size3; ii++){
        if ( (output[ii]!=complex<MyFloat> (0.,0.)) && (sum<10)  ){
            cerr<< " output alpha "<< ii<< " " << output[ii] <<endl; sum+=1;} 
    }
    //cerr <<  "output alpha" << output[0] << " " << output[1] << " " <<  output[2] << " " <<  output[3] << " " << endl;




}







#endif //IC_CONSTRAINT_HPP
