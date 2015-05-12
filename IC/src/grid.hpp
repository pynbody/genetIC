#ifndef __GRID_HPP
#define __GRID_HPP

#include <cassert>
#include "fft.hpp"

using namespace std;

template<typename T>
class Grid;

template<typename T>
class SuperSampleGrid;

template<typename T>
class Grid{
public:

    using TField = std::vector<std::complex<T>>;
    using TRealField = std::vector<T>;
    using PtrTField = std::shared_ptr<std::vector<std::complex<T>>>;

private:
    std::shared_ptr<std::vector<std::complex<T>>> pField;
    bool fieldFourier; //< is the field in k-space (true) or x-space (false)?

    // The grid offsets after Zeldovich approximation is applied
    // (nullptr before that):
    std::shared_ptr<TRealField> pOff_x;
    std::shared_ptr<TRealField> pOff_y;
    std::shared_ptr<TRealField> pOff_z;

    T hFactor, cellMass;

protected:
    T massFac;

public:

    const T simsize,boxsize,dx,x0,y0,z0;
    const size_t size;
    const size_t size2;
    const size_t size3;


    std::vector<size_t> particleArray; // just a list of particles on this grid for one purpose or another

    std::vector<T*> particleProperties; // a list of particle properties



    Grid(T simsize, size_t n, T dx=1.0, T x0=0.0, T y0=0.0, T z0=0.0) :
            dx(dx), x0(x0), y0(y0), z0(z0),
            size(n), size2(n*n), size3(n*n*n),
            boxsize(dx*n), simsize(simsize),
            massFac(1.0)
    {
        pField = std::make_shared<std::vector<std::complex<T>>>(size3,0);
        pField->shrink_to_fit();
    }

    Grid(std::true_type no_fields, T simsize, size_t n, T dx=1.0, T x0=0.0, T y0=0.0, T z0=0.0) :
            dx(dx), x0(x0), y0(y0), z0(z0),
            size(n), size2(n*n), size3(n*n*n),
            boxsize(dx*n), simsize(simsize),
            massFac(1.0)
    {

    }


    Grid(size_t n): size(n), size2(n*n), size3(n*n*n), dx(1.0), x0(0),y0(0),z0(0), boxsize(n), simsize(0), massFac(1.0) {

    }

    ///////////////////////////////////////
    //  Field manipulation routines
    ///////////////////////////////////////

    virtual TField & get_field_fourier() {
        if(!fieldFourier) {
            fft(pField->data(),pField->data(),size,1);
            fieldFourier=true;
        }
        return *pField;
    }

    virtual TField & get_field_real() {
        if(fieldFourier) {
            fft(pField->data(),pField->data(),size,-1);
            fieldFourier=false;
        }
        return *pField;
    }

    virtual TField & get_field() {
        return *pField;
    }

    auto get_offset_fields() {
        return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    virtual bool is_field_fourier()  const {
        return fieldFourier;
    }

    virtual void get_particle_no_offset(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {

        x = (*pOff_x)[id];
        y = (*pOff_y)[id];
        z = (*pOff_z)[id];

        vx = (*pOff_x)[id]*hFactor;
        vy = (*pOff_y)[id]*hFactor;
        vz = (*pOff_z)[id]*hFactor;

        cellMassi = cellMass*massFac;
        eps = dx*0.007143; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    virtual void get_particle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        get_particle_no_offset(id, x, y, z, vx, vy, vz, cellMassi, eps);
        add_centroid_location(id,x,y,z);
    }

    virtual std::shared_ptr<Grid<T>> massSplit(T massRatio) {
        massFac = massRatio;
        auto gas = std::make_shared<Grid<T>>(*this);
        gas->massFac=1.0-massRatio;
        return gas;
    }


    virtual void zeldovich(T hfac, T particlecellMass) {

        hFactor = hfac;
        cellMass = particlecellMass;

        cout << "Applying Zeldovich approximation; grid cell size=" << dx << " Mpc/h...";
        cout.flush();

        // make three arrays for manipulating in fourier space
        auto psift1k = std::vector<std::complex<T>>(size3,0);
        auto psift2k = std::vector<std::complex<T>>(size3,0);
        auto psift3k = std::vector<std::complex<T>>(size3,0);

        // get a reference to the density field in fourier space
        auto & pField_k = get_field_fourier();

        int iix, iiy, iiz;
        T kfft;
        size_t idx;

        T kw = 2.*M_PI/boxsize;

        #pragma omp parallel for schedule(static) default(shared) private(iix, iiy, iiz, kfft, idx)
        for(int ix=0; ix<size;ix++){
            for(int iy=0;iy<size;iy++){
                for(int iz=0;iz<size;iz++){

                    idx = static_cast<size_t>((ix*size+iy)*(size)+iz);

                    if( ix>size/2 ) iix = ix - size; else iix = ix;
                    if( iy>size/2 ) iiy = iy - size; else iiy = iy;
                    if( iz>size/2 ) iiz = iz - size; else iiz = iz;

                    kfft = sqrt(iix*iix+iiy*iiy+iiz*iiz);

                    psift1k[idx].real(-pField_k[idx].imag()/(T)(kfft*kfft)*iix/kw);
                    psift1k[idx].imag(pField_k[idx].real()/(T)(kfft*kfft)*iix/kw);
                    psift2k[idx].real(-pField_k[idx].imag()/(T)(kfft*kfft)*iiy/kw);
                    psift2k[idx].imag(pField_k[idx].real()/(T)(kfft*kfft)*iiy/kw);
                    psift3k[idx].real(-pField_k[idx].imag()/(T)(kfft*kfft)*iiz/kw);
                    psift3k[idx].imag(pField_k[idx].real()/(T)(kfft*kfft)*iiz/kw);
                }
            }
        }

        psift1k[0]=complex<T>(0.,0.);
        psift2k[0]=complex<T>(0.,0.);
        psift3k[0]=complex<T>(0.,0.);

        fft(psift1k.data(),psift1k.data(),size,-1); //the output .imag() part is non-zero because of the Nyquist frequency, but this is not used anywhere else
        fft(psift2k.data(),psift2k.data(),size,-1); //same
        fft(psift3k.data(),psift3k.data(),size,-1); //same

        pOff_x = std::make_shared<std::vector<T>>(size3,0);
        pOff_y = std::make_shared<std::vector<T>>(size3,0);
        pOff_z = std::make_shared<std::vector<T>>(size3,0);


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


    ///////////////////////////////////////
    //  Index manipulation routines
    ///////////////////////////////////////

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

    T get_abs_k_coordinates(long id) const {
        int x,y,z;
        get_k_coordinates(id,x,y,z);
        return sqrt(x*x+y*y+z*z);
    }


    void get_centroid_location(size_t id, T &xc, T &yc, T &zc) const {
        int x, y, z;
        get_coordinates(id,x,y,z);
        xc = x0+x*dx+dx/2;
        yc = y0+y*dx+dx/2;
        zc = z0+z*dx+dx/2;

    }

    void add_centroid_location(size_t id, T &xc, T &yc, T &zc) const {
        int x, y, z;
        get_coordinates(id,x,y,z);
        xc += x0+x*dx+dx/2;
        yc += y0+y*dx+dx/2;
        zc += z0+z*dx+dx/2;

        // always wrap at the BASE level:
        x = fmod(x,simsize);
        if(x<0) x+=simsize;
        y = fmod(y,simsize);
        if(y<0) y+=simsize;
        z = fmod(z,simsize);
        if(z<0) z+=simsize;

    }

    tuple<T, T, T> get_centroid_location(long id) const {
        T xc,yc,zc;
        get_centroid_location(id,xc,yc,zc);
        return std::make_tuple(xc,yc,zc);
    }

    vector<size_t> get_ids_in_cube(T x0c, T y0c, T z0c, T dxc) {
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




};


template<typename T>
class SuperSampleGrid : public Grid<T> {
private:
    std::shared_ptr<Grid<T>> pUnderlying;
    int factor;
    int factor3;

protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;

public:
    SuperSampleGrid(std::shared_ptr<Grid<T>> pUnderlying, int factor):
            Grid<T>(
                    pUnderlying->simsize, pUnderlying->size*factor,
                    pUnderlying->dx, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0)
    {

        this->massFac = 1.0;
        factor3=factor*factor*factor;
    }

    // all the field manipulation routines MUST NOT be called, since we are
    // going to interpolate on the fly
    virtual TField & get_field_fourier() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    virtual TField & get_field_real() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    virtual TField & get_field() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    /*
    virtual std::tuple<TRealField &, TRealField &, TRealField &> get_offset_fields() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }
    */

    virtual bool is_field_fourier()  const override  {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }


    virtual void get_particle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        size_t id_underlying = id/factor3;

        pUnderlying->get_particle_no_offset(id_underlying, x, y, z, vx, vy, vz, cellMassi, eps);
        this->add_centroid_location(id,x,y,z);

        cellMassi*=this->massFac/factor3;

    }

    virtual std::shared_ptr<Grid<T>> massSplit(T massRatio) override {
        cerr << "WARNING: massSplit has been called on a supersampled grid. This is unlikely to be what you want...?" << endl;
        this->massFac = massRatio;
        auto gas = std::make_shared<SuperSampleGrid<T>>(this->pUnderlying, factor);
        gas->massFac=1.0-massRatio;
        return gas;
    }

    virtual void zeldovich(T hfac, T particlecellMass) override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

};




#endif
