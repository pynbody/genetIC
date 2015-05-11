#ifndef __GRID_HPP
#define __GRID_HPP

#include <cassert>
#include "fft.hpp"

using namespace std;

template<typename MyFloat>
struct grid_struct{

    long grid[3];
    long coords[3];
    MyFloat absval;
    MyFloat delta;

};

template<typename MyFloat>
class Grid{
private:
    std::shared_ptr<std::vector<std::complex<MyFloat>>> pField;
    bool fieldFourier; //< is the field in k-space (true) or x-space (false)?

    // The grid offsets after Zeldovich approximation is applied
    // (nullptr before that):
    std::shared_ptr<std::vector<MyFloat>> pOff_x;
    std::shared_ptr<std::vector<MyFloat>> pOff_y;
    std::shared_ptr<std::vector<MyFloat>> pOff_z;

    MyFloat hFactor, cellMass;

protected:
    MyFloat massFac;

public:

    const MyFloat simsize,boxsize,dx,x0,y0,z0;
    const size_t size;
    const size_t size2;
    const size_t size3;

    std::vector<size_t> particleArray; // just a list of particles on this grid for one purpose or another

    std::vector<MyFloat*> particleProperties; // a list of particle properties



    Grid(MyFloat simsize, size_t n, MyFloat dx=1.0, MyFloat x0=0.0, MyFloat y0=0.0, MyFloat z0=0.0) :
            dx(dx), x0(x0), y0(y0), z0(z0),
            size(n), size2(n*n), size3(n*n*n),
            boxsize(dx*n), simsize(simsize),
            massFac(1.0)
    {
        pField = std::make_shared<std::vector<std::complex<MyFloat>>>(size3,0);
        pField->shrink_to_fit();
    }

    Grid(size_t n): size(n), size2(n*n), size3(n*n*n), dx(1.0), x0(0),y0(0),z0(0), boxsize(n), simsize(0), massFac(1.0) {

    }

    std::vector<std::complex<MyFloat>> & get_field_fourier() {
        if(!fieldFourier) {
            fft(pField->data(),pField->data(),size,1);
            fieldFourier=true;
        }
        return *pField;
    }

    std::vector<std::complex<MyFloat>> & get_field_real() {
        if(fieldFourier) {
            fft(pField->data(),pField->data(),size,-1);
            fieldFourier=false;
        }
        return *pField;
    }

    std::vector<std::complex<MyFloat>> & get_field() {
        return *pField;
    }

    auto get_offset_fields() {
        return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    bool is_field_fourier()  const {
        return fieldFourier;
    }


    long find_next_ind(long index, const int step[3]) const
    {
        int grid[3];
        std::tie(grid[0],grid[1],grid[2])=get_coordinates(index);

        grid[0]+=step[0];
        grid[1]+=step[1];
        grid[2]+=step[2];

        return this->get_index(grid); // N.B. does wrapping inside get_index
    }

    long find_next_ind_no_wrap(long index, const int step[3]) const
    {

        int grid[3];
        std::tie(grid[0],grid[1],grid[2])=get_coordinates(index);

        grid[0]+=step[0];
        grid[1]+=step[1];
        grid[2]+=step[2];


        return this->get_index_no_wrap(grid);
    }

    void wrap(int &x, int &y, int &z) const
    {
        x = x%size;
        y = y%size;
        z = z%size;
        if(x<0) x+=size;
        if(y<0) y+=size;
        if(z<0) z+=size;
    }

    void wrap(int pos[3]) const
    {
        wrap(pos[0],pos[1],pos[2]);
    }


    long get_index(int x, int y, int z) const
    {

        long size=this->size;

        wrap(x,y,z);

        long index=(x*size+y);
        index*=size;
        index+=z;

        return index;

    }

    long get_index_no_wrap(int x, int y, int z) const
    {

        long size=this->size;

        if(x<0 || x>=size || y<0 || y>=size || z<0 || z>=size)
            throw std::runtime_error("Grid index out of range in get_index_no_wrap");

        long index=(x*size+y);
        index*=size;
        index+=z;

        return index;

    }

    long get_index(const int pos[3]) const {
        return get_index(pos[0], pos[1], pos[2]);
    }

    long get_index_no_wrap(const int pos[3]) const {
        return get_index_no_wrap(pos[0], pos[1], pos[2]);
    }

    void get_coordinates(long id, int &x, int &y, int &z) const {
        x = (int) (id/size2);
        y = (int) (id%size2)/size;
        z = (int) (id%size);

        // TODO: optimization - following check should be removed at some point:
        if(get_index(x,y,z)!=id) {
            cerr << "ERROR in get_coordinates";
            cerr << "id=" << id << " x,y,z=" << x << "," << y << "," << z << endl;
            cerr << "which gives " << get_index(x,y,z) << endl;
            assert(false);
        }
    }

    tuple<int, int, int> get_coordinates(long id) const {
        int x, y, z;

        get_coordinates(id, x,y,z);

        return std::make_tuple(x,y,z);
    }

    void get_k_coordinates(long id, int &x, int &y, int &z) const {
        get_coordinates(id,x,y,z);
        if(x>size/2) x=x-size;
        if(y>size/2) y=y-size;
        if(z>size/2) z=z-size;
    }

    tuple<int, int, int> get_k_coordinates(long id) const {
        int x, y, z;

        get_k_coordinates(id, x,y,z);

        return std::make_tuple(x,y,z);
    }

    MyFloat get_abs_k_coordinates(long id) const {
        int x,y,z;
        get_k_coordinates(id,x,y,z);
        return sqrt(x*x+y*y+z*z);
    }


    void get_centroid_location(size_t id, MyFloat &xc, MyFloat &yc, MyFloat &zc) const {
        int x, y, z;
        get_coordinates(id,x,y,z);
        xc = x0+x*dx+dx/2;
        yc = y0+y*dx+dx/2;
        zc = z0+z*dx+dx/2;

    }

    tuple<MyFloat, MyFloat, MyFloat> get_centroid_location(long id) const {
        MyFloat xc,yc,zc;
        get_centroid_location(id,xc,yc,zc);
        return std::make_tuple(xc,yc,zc);
    }

    vector<size_t> get_ids_in_cube(MyFloat x0c, MyFloat y0c, MyFloat z0c, MyFloat dxc) {
        // return all the grid IDs whose centres lie within the specified cube
        vector<size_t> ids;

        // TODO: optimization, set the storage size of ids here.

        int xa=((int) floor((x0c-x0-dxc/2+dx/2)/dx));
        int ya=((int) floor((y0c-y0-dxc/2+dx/2)/dx));
        int za=((int) floor((z0c-z0-dxc/2+dx/2)/dx));

        int xb=((int) floor((x0c-x0+dxc/2-dx/2)/dx));
        int yb=((int) floor((y0c-y0+dxc/2-dx/2)/dx));
        int zb=((int) floor((z0c-z0+dxc/2-dx/2)/dx));

        for(int x=xa; x<=xb; x++) {
            for(int y=ya; y<=yb; y++) {
                for(int z=za; z<=zb; z++) {
                    ids.push_back(get_index(x,y,z));
                }
            }
        }

        return ids;

    }

    void get_particle(size_t id, MyFloat &x, MyFloat &y, MyFloat &z, MyFloat &vx, MyFloat &vy, MyFloat &vz, MyFloat &cellMassi, MyFloat &eps) const
    {
        get_centroid_location(id,x,y,z);
        x += (*pOff_x)[id];
        y += (*pOff_y)[id];
        z += (*pOff_z)[id];

        // always wrap at the BASE level:
        x = fmod(x,simsize);
        if(x<0) x+=simsize;
        y = fmod(y,simsize);
        if(y<0) y+=simsize;
        z = fmod(z,simsize);
        if(z<0) z+=simsize;



        vx = (*pOff_x)[id]*hFactor;
        vy = (*pOff_y)[id]*hFactor;
        vz = (*pOff_z)[id]*hFactor;

        cellMassi = cellMass*massFac;
        eps = dx*0.007143; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    auto massSplit(MyFloat massRatio) {
        massFac = massRatio;
        auto gas = std::make_shared<Grid<MyFloat>>(*this);
        gas->massFac=1.0-massRatio;
        return gas;
    }

    void zeldovich(MyFloat hfac, MyFloat particlecellMass) {

        hFactor = hfac;
        cellMass = particlecellMass;

        cout << "Applying Zeldovich approximation; grid cell size=" << dx << " Mpc/h...";
        cout.flush();

        // make three arrays for manipulating in fourier space
        auto psift1k = std::vector<std::complex<MyFloat>>(size3,0);
        auto psift2k = std::vector<std::complex<MyFloat>>(size3,0);
        auto psift3k = std::vector<std::complex<MyFloat>>(size3,0);

        // get a reference to the density field in fourier space
        auto & pField_k = get_field_fourier();

        int iix, iiy, iiz;
        MyFloat kfft;
        size_t idx;

        MyFloat kw = 2.*M_PI/boxsize;

        #pragma omp parallel for schedule(static) default(shared) private(iix, iiy, iiz, kfft, idx)
        for(int ix=0; ix<size;ix++){
            for(int iy=0;iy<size;iy++){
                for(int iz=0;iz<size;iz++){

                    idx = static_cast<size_t>((ix*size+iy)*(size)+iz);

                    if( ix>size/2 ) iix = ix - size; else iix = ix;
                    if( iy>size/2 ) iiy = iy - size; else iiy = iy;
                    if( iz>size/2 ) iiz = iz - size; else iiz = iz;

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

        psift1k[0]=complex<MyFloat>(0.,0.);
        psift2k[0]=complex<MyFloat>(0.,0.);
        psift3k[0]=complex<MyFloat>(0.,0.);

        fft(psift1k.data(),psift1k.data(),size,-1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
        fft(psift2k.data(),psift2k.data(),size,-1); //same
        fft(psift3k.data(),psift3k.data(),size,-1); //same

        pOff_x = std::make_shared<std::vector<MyFloat>>(size3,0);
        pOff_y = std::make_shared<std::vector<MyFloat>>(size3,0);
        pOff_z = std::make_shared<std::vector<MyFloat>>(size3,0);


        //apply ZA:
        #pragma omp parallel for schedule(static) default(shared) private(iix, iiy, iiz, kfft, idx)
        for(int ix=0;ix<size;ix++) {
            for(int iy=0;iy<size;iy++) {
                for(int iz=0;iz<size;iz++) {

                    idx = static_cast<size_t>((ix*size+iy)*(size)+iz);

                    // position offset in physical coordinates
                    (*pOff_x)[idx] = psift1k[idx].real();
                    (*pOff_y)[idx] = psift2k[idx].real();
                    (*pOff_z)[idx] = psift3k[idx].real();

                }
            }
        }

        cout << "done."<<endl;

    }

/*
    void add_grid(MyFloat *Pos1, MyFloat *Pos2, MyFloat *Pos3, MyFloat boxlen) {
        if(boxlen<0)
            boxlen = dx*size;

        MyFloat Mean1=0, Mean2=0, Mean3=0;
        size_t idx;

        for(int ix=0;ix<size;ix++) {
            for(int iy=0;iy<size;iy++) {
                for(int iz=0;iz<size;iz++) {

                    idx = static_cast<size_t>((ix*size+iy)*size+iz);

                    // position in physical coordinates
                    Pos1[idx]+= ix*dx+dx/2+x0;
                    Pos2[idx]+= iy*dx+dx/2+y0;
                    Pos3[idx]+= iz*dx+dx/2+z0;

                    // always wrap at the BASE level:
                    Pos1[idx] = fmod(Pos1[idx],boxlen);
                    if(Pos1[idx]<0) Pos1[idx]+=boxlen;
                    Pos2[idx] = fmod(Pos2[idx],boxlen);
                    if(Pos2[idx]<0) Pos2[idx]+=boxlen;
                    Pos3[idx] = fmod(Pos3[idx],boxlen);
                    if(Pos3[idx]<0) Pos3[idx]+=boxlen;


                    Mean1+=Pos1[idx];
                    Mean2+=Pos2[idx];
                    Mean3+=Pos3[idx];

                }
            }
        }


        cout<< "Box/2="<< boxlen/2.<< " Mpc/h, Mean position x,y,z: "<< Mean1/(MyFloat(size*size*size))<<" "<< Mean2/(MyFloat(size*size*size))<<" "<<Mean3/(MyFloat(size*size*size))<< " Mpc/h"<<  endl;

    }
    */


};

#endif
