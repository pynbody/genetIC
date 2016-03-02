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
        grid.getCentroidLocation(index, xp, yp, zp);

        xp=grid.getWrappedDelta(xp,x0);
        yp=grid.getWrappedDelta(yp,y0);
        zp=grid.getWrappedDelta(zp,z0);

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

        else{cerr<< "Wrong value for parameter 'direc' in function 'cen_deriv4_alpha'."<< endl; exit(1);}

        for(int di=0; di<3; di++) {
            long ind_p1, ind_m1, ind_p2, ind_m2;
            int step1[3]={0,0,0};
            int neg_step1[3]={0,0,0};
            step1[di]=1;
            neg_step1[di]=-1;

            // N.B. can't wrap - might be on subgrid -> we do it anyway! this will likely break with zoom-in TODO
            //ind_m1=grid.findNextIndNoWrap(index, neg_step1);
            //ind_p1=grid.findNextIndNoWrap(index, step1);
            //ind_m2=grid.findNextIndNoWrap(ind_m1, neg_step1);
            //ind_p2=grid.findNextIndNoWrap(ind_p1, step1);

            ind_m1=grid.findNextInd(index, neg_step1);
            ind_p1=grid.findNextInd(index, step1);
            ind_m2=grid.findNextInd(ind_m1, neg_step1);
            ind_p2=grid.findNextInd(ind_p1, step1);

            MyFloat a=-1./12./grid.dx, b=2./3./grid.dx;  //the signs here so that L ~ - Nabla Phi

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

        for(size_t i=0;i<particleArray.size();i++) {
            grid.getCentroidLocation(particleArray[i],xb,yb,zb);

            x0+=grid.getWrappedDelta(xb,xa);
            y0+=grid.getWrappedDelta(yb,ya);
            z0+=grid.getWrappedDelta(zb,za);
        }

        x0/=particleArray.size();
        y0/=particleArray.size();
        z0/=particleArray.size();
        x0+=xa;
        y0+=ya;
        z0+=za;

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

//OLD angular momentum class, temporary fix below (because I couldn't figure out how to pass the direction):
    void angmom(int direction) {
        throw std::runtime_error("Angmom cannot run - x0,y0,z0 coordinates not yet passed to new constraint class");

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

    void angmom0() { //L_x

        MyFloat x0, y0, z0;
        x0=this->x0;
        y0=this->y0;
        z0=this->z0;
  
        for(size_t i=0;i<particleArray.size();i++) {
            cen_deriv4_alpha(particleArray[i], 0, x0, y0, z0);
        }

        fft(output.data(), output.data(), grid.size, 1);
        // The constraint as derived is on the potential. By considering
        // unitarity of FT, we can FT the constraint to get the constraint
        // on the density.
        poiss(output.data(), output.data(), grid.size, grid.boxsize,
              cosmology.scalefactor, cosmology.OmegaM0);

    }

    void angmom1() { //L_y
        
        MyFloat x0, y0, z0;
        x0=this->x0;
        y0=this->y0;
        z0=this->z0;

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

    void angmom2() { //L_z
        
        MyFloat x0, y0, z0;
        x0=this->x0;
        y0=this->y0;
        z0=this->z0;

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
    //temporary fix to calculate angular momentum (by Nina): (how to pass the direction here?)
    dispatch.add_class_route("L0",&CC::angmom0);
    dispatch.add_class_route("L1",&CC::angmom1);
    dispatch.add_class_route("L2",&CC::angmom2);

    dispatch.run(name);

}







#endif //IC_CONSTRAINT_HPP
