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
class SubSampleGrid;

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



    Grid(T simsize, size_t n, T dx=1.0, T x0=0.0, T y0=0.0, T z0=0.0, bool withField=true) :
            massFac(1.0),
            simsize(simsize), boxsize(dx*n),
            dx(dx), x0(x0), y0(y0), z0(z0),
            size(n), size2(n*n), size3(n*n*n)
    {
        if(withField) {
            pField = std::make_shared<std::vector<std::complex<T>>>(size3,0);
            pField->shrink_to_fit();
        }
    }




    Grid(size_t n): massFac(1.0),  simsize(0), boxsize(n),
                    dx(1.0), x0(0),y0(0),z0(0), size(n), size2(n*n), size3(n*n*n)
    {

    }

    ///////////////////////////////////////
    //  Field manipulation routines
    ///////////////////////////////////////

    virtual TField & getFieldFourier() {
        assert(pField!=nullptr);
        if(!fieldFourier) {
            fft(pField->data(),pField->data(),size,1);
            fieldFourier=true;
        }
        return *pField;
    }

    virtual TField & getFieldReal() {
        assert(pField!=nullptr);
        if(fieldFourier) {
            fft(pField->data(),pField->data(),size,-1);
            fieldFourier=false;
        }
        return *pField;
    }

    virtual TField & getField() {
        assert(pField!=nullptr);
        return *pField;
    }

    auto getOffsetFields() {
        return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    virtual bool isFieldFourier()  const {
        return fieldFourier;
    }

    virtual void getParticle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        getParticleNoOffset(id, x, y, z, vx, vy, vz, cellMassi, eps);
        addCentroidLocation(id,x,y,z);
    }

    virtual T getMass() const {
        return cellMass*massFac;
    }

    virtual T getEps() const {
        return dx*0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    virtual void getParticleNoOffset(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {

        x = (*pOff_x)[id];
        y = (*pOff_y)[id];
        z = (*pOff_z)[id];

        vx = (*pOff_x)[id]*hFactor;
        vy = (*pOff_y)[id]*hFactor;
        vz = (*pOff_z)[id]*hFactor;

        cellMassi = getMass();
        eps = getEps();
    }

    virtual void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {

        vx=getFieldInterpolated(x,y,z,*pOff_x);
        vy=getFieldInterpolated(x,y,z,*pOff_y);
        vz=getFieldInterpolated(x,y,z,*pOff_z);
        x+=vx;
        y+=vy;
        z+=vz;

        simWrap(x,y,z);

        vx*=hFactor;
        vy*=hFactor;
        vz*=hFactor;

        cellMassi = getMass();
        eps = getEps();
    }




    template<typename TField>
    typename TField::value_type getFieldInterpolated(const T &x, const T &y, const T &z, const TField & pField) const {

        int x_p_0, y_p_0, z_p_0, x_p_1, y_p_1, z_p_1;

        // grid coordinates of parent cell starting to bottom-left
        // of our current point
        x_p_0 = (int) floor(((x-x0)/dx - 0.5));
        y_p_0 = (int) floor(((y-y0)/dx - 0.5));
        z_p_0 = (int) floor(((z-z0)/dx - 0.5));

        // grid coordinates of top-right
        x_p_1 = x_p_0+1;
        y_p_1 = y_p_0+1;
        z_p_1 = z_p_0+1;

        // weights, which are the distance to the centre point of the
        // upper-right cell, in grid units (-> maximum 1)

        T xw0,yw0,zw0,xw1,yw1,zw1;
        xw0 = ((T) x_p_1 + 0.5) - ((x-x0)/dx);
        yw0 = ((T) y_p_1 + 0.5) - ((y-y0)/dx);
        zw0 = ((T) z_p_1 + 0.5) - ((z-z0)/dx);

        xw1 = 1.-xw0;
        yw1 = 1.-yw0;
        zw1 = 1.-zw0;

        assert(xw0<=1.0 && xw0>=0.0);

        // allow things on the boundary to 'saturate' value, but beyond boundary
        // is not acceptable
        //
        // TODO - in some circumstances we may wish to replace this with wrapping
        // but not all circumstances!
        int size_i = static_cast<int>(size);
        assert(x_p_1<=size_i);
        if(x_p_1==size_i) x_p_1=size_i-1;
        assert(y_p_1<=size_i);
        if(y_p_1==size_i) y_p_1=size_i-1;
        assert(z_p_1<=size_i);
        if(z_p_1==size_i) z_p_1=size_i-1;

        assert(x_p_0>=-1);
        if(x_p_0==-1) x_p_0=0;
        assert(y_p_0>=-1);
        if(y_p_0==-1) y_p_0=0;
        assert(z_p_0>=-1);
        if(z_p_0==-1) z_p_0=0;


        return xw0*yw0*zw1*pField[getIndexNoWrap(x_p_0,y_p_0,z_p_1)] +
               xw1*yw0*zw1*pField[getIndexNoWrap(x_p_1,y_p_0,z_p_1)] +
               xw0*yw1*zw1*pField[getIndexNoWrap(x_p_0,y_p_1,z_p_1)] +
               xw1*yw1*zw1*pField[getIndexNoWrap(x_p_1,y_p_1,z_p_1)] +
               xw0*yw0*zw0*pField[getIndexNoWrap(x_p_0,y_p_0,z_p_0)] +
               xw1*yw0*zw0*pField[getIndexNoWrap(x_p_1,y_p_0,z_p_0)] +
               xw0*yw1*zw0*pField[getIndexNoWrap(x_p_0,y_p_1,z_p_0)] +
               xw1*yw1*zw0*pField[getIndexNoWrap(x_p_1,y_p_1,z_p_0)] ;
    }




    virtual std::shared_ptr<Grid<T>> massSplit(T massRatio) {
        massFac = 1.0-massRatio;
        auto gas = std::make_shared<Grid<T>>(*this);
        gas->massFac=massRatio;
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
        auto & pField_k = getFieldFourier();

        int iix, iiy, iiz;
        T kfft;
        size_t idx;

        T kw = 2.*M_PI/boxsize;

        #pragma omp parallel for schedule(static) default(shared) private(iix, iiy, iiz, kfft, idx)
        for(size_t ix=0; ix<size;ix++){
            for(size_t iy=0;iy<size;iy++){
                for(size_t iz=0;iz<size;iz++){

                    idx = (ix*size+iy)*size+iz;

                    if( ix>size/2 ) iix = static_cast<int>(ix) - size; else iix = ix;
                    if( iy>size/2 ) iiy = static_cast<int>(iy) - size; else iiy = iy;
                    if( iz>size/2 ) iiz = static_cast<int>(iz) - size; else iiz = iz;

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
        for(size_t ix=0;ix<size;ix++) {
            for(size_t iy=0;iy<size;iy++) {
                for(size_t iz=0;iz<size;iz++) {

                    idx = (ix*size+iy)*size+iz;

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

    long findNextInd(long index, const int step[3]) const
    {
        int grid[3];
        std::tie(grid[0],grid[1],grid[2])=getCoordinates(index);

        grid[0]+=step[0];
        grid[1]+=step[1];
        grid[2]+=step[2];

        return this->getIndex(grid); // N.B. does wrapping inside getIndex
    }

    long findNextIndNoWrap(long index, const int step[3]) const
    {

        int grid[3];
        std::tie(grid[0],grid[1],grid[2])=getCoordinates(index);

        grid[0]+=step[0];
        grid[1]+=step[1];
        grid[2]+=step[2];


        return this->getIndexNoWrap(grid);
    }

    void wrap(int &x, int &y, int &z) const
    {
#ifdef SAFER_SLOWER
        x = x%size;
        y = y%size;
        z = z%size;
#else
        if(x>size) x-=size;
        if(y>size) y-=size;
        if(z>size) z-=size;
#endif
        if(x<0) x+=size;
        if(y<0) y+=size;
        if(z<0) z+=size;
    }

    void simWrap(T &x, T &y, T &z) const
    {
        x = fmod(x,simsize);
        if(x<0) x+=simsize;
        y = fmod(y,simsize);
        if(y<0) y+=simsize;
        z = fmod(z,simsize);
        if(z<0) z+=simsize;
    }

    void wrap(int pos[3]) const
    {
        wrap(pos[0],pos[1],pos[2]);
    }


    long getIndex(int x, int y, int z) const
    {

        long size=this->size;

        wrap(x,y,z);

        long index=(x*size+y);
        index*=size;
        index+=z;

        return index;

    }

    size_t getIndexNoWrap(size_t x, size_t y, size_t z) const
    {

#ifdef SAFER_SLOWER
        if(x<0 || x>=size || y<0 || y>=size || z<0 || z>=size)
            throw std::runtime_error("Grid index out of range in getIndexNoWrap");
#endif
        return (x*size+y)*size+z;
    }

    long getIndex(const int pos[3]) const {
        return getIndex(pos[0], pos[1], pos[2]);
    }

    long getIndexNoWrap(const int pos[3]) const {
        return getIndexNoWrap(pos[0], pos[1], pos[2]);
    }

    void getCoordinates(long id, int &x, int &y, int &z) const {
        x = (int) (id/size2);
        y = (int) (id%size2)/size;
        z = (int) (id%size);

        // TODO: optimization - following check should be removed at some point:
        if(getIndex(x,y,z)!=id) {
            cerr << "ERROR in getCoordinates";
            cerr << "id=" << id << " x,y,z=" << x << "," << y << "," << z << endl;
            cerr << "which gives " << getIndex(x,y,z) << endl;
            assert(false);
        }
    }

    tuple<int, int, int> getCoordinates(long id) const {
        int x, y, z;

        getCoordinates(id, x,y,z);

        return std::make_tuple(x,y,z);
    }

    void getKCoordinates(long id, int &x, int &y, int &z) const {
        getCoordinates(id,x,y,z);
        if(x>size/2) x=x-size;
        if(y>size/2) y=y-size;
        if(z>size/2) z=z-size;
    }

    tuple<int, int, int> getKCoordinates(long id) const {
        int x, y, z;

        getKCoordinates(id, x,y,z);

        return std::make_tuple(x,y,z);
    }

    T getAbsKCoordinates(long id) const {
        int x,y,z;
        getKCoordinates(id,x,y,z);
        return sqrt(x*x+y*y+z*z);
    }


    void getCentroidLocation(size_t id, T &xc, T &yc, T &zc) const {
        int x, y, z;
        getCoordinates(id,x,y,z);
        xc = x0+x*dx+dx/2;
        yc = y0+y*dx+dx/2;
        zc = z0+z*dx+dx/2;

    }

    void addCentroidLocation(size_t id, T &xc, T &yc, T &zc) const {
        int x, y, z;
        getCoordinates(id,x,y,z);
        xc += x0+x*dx+dx/2;
        yc += y0+y*dx+dx/2;
        zc += z0+z*dx+dx/2;

        // always wrap at the BASE level:
        simWrap(xc,yc,zc);

    }

    tuple<T, T, T> getCentroidLocation(long id) const {
        T xc,yc,zc;
        getCentroidLocation(id,xc,yc,zc);
        return std::make_tuple(xc,yc,zc);
    }

    vector<size_t> getIdsInCube(T x0c, T y0c, T z0c, T dxc) {
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
                    ids.push_back(getIndex(x,y,z));
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
                    pUnderlying->dx/factor, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0, false),
            pUnderlying(pUnderlying), factor(factor)
    {

        this->massFac = 1.0;
        factor3=factor*factor*factor;
    }

    // all the field manipulation routines MUST NOT be called, since we are
    // going to interpolate on the fly
    virtual TField & getFieldFourier() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    virtual TField & getFieldReal() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    virtual TField & getField() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    /*
    virtual std::tuple<TRealField &, TRealField &, TRealField &> getOffsetFields() override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }
    */

    virtual bool isFieldFourier()  const override  {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

    virtual T getMass() const override {
        return pUnderlying->getMass()*this->massFac/factor3;
    }


    virtual void getParticle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        // get location
        this->getCentroidLocation(id,x,y,z);

        // interpolate from underlying grid
        pUnderlying->getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);

        // adjust mass
        cellMassi*=this->massFac/factor3;

    }

    virtual std::shared_ptr<Grid<T>> massSplit(T massRatio) override {
        cerr << "WARNING: massSplit has been called on a supersampled grid. This is unlikely to be what you want...?" << endl;
        this->massFac = 1.0-massRatio;
        auto gas = std::make_shared<SuperSampleGrid<T>>(this->pUnderlying, factor);
        gas->massFac=massRatio;
        return gas;
    }

    virtual void zeldovich(T hfac, T particlecellMass) override {
        throw std::runtime_error("SuperSampleGrid - does not contain an actual field in memory");
    }

};



template<typename T>
class SubSampleGrid : public Grid<T> {
private:
    std::shared_ptr<Grid<T>> pUnderlying;
    int factor;
    int factor3;

protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;

public:
    SubSampleGrid(std::shared_ptr<Grid<T>> pUnderlying, int factor):
            Grid<T>(
                    pUnderlying->simsize, pUnderlying->size/factor,
                    pUnderlying->dx*factor, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0, false),
            pUnderlying(pUnderlying), factor(factor)
    {
        if(pUnderlying->size%factor!=0)
          throw std::runtime_error("SubSampleGrid - factor must be a divisor of the original grid size");

        this->massFac = 1.0;
        factor3=factor*factor*factor;
    }

    // all the field manipulation routines MUST NOT be called, since we are
    // going to interpolate on the fly
    virtual TField & getFieldFourier() override {
        throw std::runtime_error("SubSampleGrid - does not contain an actual field in memory");
    }

    virtual TField & getFieldReal() override {
        throw std::runtime_error("SubSampleGrid - does not contain an actual field in memory");
    }

    virtual TField & getField() override {
        throw std::runtime_error("SubSampleGrid - does not contain an actual field in memory");
    }

    /*
    virtual std::tuple<TRealField &, TRealField &, TRealField &> getOffsetFields() override {
        throw std::runtime_error("SubSampleGrid - does not contain an actual field in memory");
    }
    */

    virtual bool isFieldFourier()  const override  {
        throw std::runtime_error("SubSampleGrid - does not contain an actual field in memory");
    }

    virtual T getMass() const override {
        return pUnderlying->getMass()*this->massFac/factor3;
    }


    virtual void getParticle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
      int x0,y0,z0;
      this->getCoordinates(id, x0, y0, z0);
      auto x1=x0+factor, y1=y0+factor, z1=z0+factor;

      T xt,yt,zt,vxt,vyt,vzt,cellMassit,epst;

      x=0;
      y=0;
      z=0;
      vx=0;
      vy=0;
      vz=0;
      cellMassi=0;
      eps=0;

      // construct our virtual values from average over the underlying cells
      for(auto xi=x0; xi<x1; ++xi) {
        for (auto yi=y0; yi<y1; ++yi) {
          for (auto zi=z0; zi<z1; ++zi) {
            pUnderlying->getParticle(pUnderlying->getIndexNoWrap(xi,yi,zi),
                                     xt,yt,zt,vxt,vyt,vzt,cellMassit,epst);
            x+=xt/factor3;
            y+=yt/factor3;
            z+=zt/factor3;
            vx+=vxt/factor3;
            vy+=vyt/factor3;
            vz+=vzt/factor3;
            eps+=epst/factor3;

            // accumulate, don't average the mass!
            cellMassi+=cellMassit;
          }
        }
      }

    }

    virtual std::shared_ptr<Grid<T>> massSplit(T massRatio) override {
        cerr << "WARNING: massSplit has been called on a subsampled grid. This is unlikely to be what you want...?" << endl;
        this->massFac = 1.0-massRatio;
        auto gas = std::make_shared<SubSampleGrid<T>>(this->pUnderlying, factor);
        gas->massFac=massRatio;
        return gas;
    }

    virtual void zeldovich(T hfac, T particlecellMass) override {
        throw std::runtime_error("SubSampleGrid - does not contain an actual field in memory");
    }

};



#endif
