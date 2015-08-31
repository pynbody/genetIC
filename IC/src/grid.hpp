#ifndef __GRID_HPP
#define __GRID_HPP

#include <cassert>
#include <set>
#include <type_traits>
#include "fft.hpp"
#include "filter.hpp"

using namespace std;

template<typename T>
class Grid;

template<typename T>
class SuperSampleGrid;

template<typename T>
class SubSampleGrid;

template<typename T>
class OffsetGrid;

template<typename T>
class SectionOfGrid;

template<typename T>
class MassScaledGrid;

template<typename T>
class Grid : public std::enable_shared_from_this<Grid<T>>  {
public:

    using TField = std::vector<std::complex<T>>;
    using TRealField = std::vector<T>;
    using PtrTField = std::shared_ptr<std::vector<std::complex<T>>>;
    using GridPtrType = std::shared_ptr<Grid<T>>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;

private:
    std::shared_ptr<std::vector<std::complex<T>>> pField;
    bool fieldFourier; //< is the field in k-space (true) or x-space (false)?

    // The grid offsets after Zeldovich approximation is applied
    // (nullptr before that):
    std::shared_ptr<TRealField> pOff_x;
    std::shared_ptr<TRealField> pOff_y;
    std::shared_ptr<TRealField> pOff_z;

    T hFactor, cellMass;
    T kMin, kMinSquared;

    std::vector<size_t> particleArray; // just a list of particles on this grid for one purpose or another
    std::vector<T*> particleProperties; // a list of particle properties

protected:

    static void upscaleParticleList(const std::vector<size_t> sourceArray,
                                    std::vector<size_t> & targetArray,
                                    const Grid<T> *source,
                                    const Grid<T> *target) {
      int x0,y0,z0,x1,y1,z1,x,y,z;

      assert(target->size>=source->size);
      assert((target->size)%(source->size)==0);
      size_t factor = target->size/source->size;
      targetArray.clear();

      for(auto id: sourceArray) {
        source->getCoordinates(id,x0,y0,z0);
        x0*=factor;
        y0*=factor;
        z0*=factor;
        x1 = x0+factor;
        y1 = y0+factor;
        z1 = z0+factor;
        for(x=x0;x<x1;++x) {
          for(y=y0;y<y1;++y) {
            for(z=z0;z<z1;++z) {
              targetArray.push_back(target->getIndexNoWrap(x,y,z));
            }
          }
        }
      }
    }

    static void downscaleParticleList(const std::vector<size_t> sourceArray,
                                    std::vector<size_t> & targetArray,
                                    const Grid<T> *source,
                                    const Grid<T> *target) {


      std::set<size_t> targetSet;
      int x,y,z;

      assert(source->size>=target->size);
      assert((source->size)%(target->size)==0);
      size_t factor = source->size/target->size;

      for(auto id: sourceArray) {
        source->getCoordinates(id,x,y,z);
        targetSet.insert(target->getIndexNoWrap(x/factor, y/factor, z/factor));
      }
      targetArray.clear();
      targetArray.insert(targetArray.end(), targetSet.begin(), targetSet.end());
    }

    void setKmin() {
        kMin = 2.*M_PI/boxsize;
        kMinSquared = kMin*kMin;
    }

public:

    const T simsize,boxsize,dx,x0,y0,z0;
    const size_t size,size2,size3;






    Grid(T simsize, size_t n, T dx=1.0, T x0=0.0, T y0=0.0, T z0=0.0, bool withField=true) :
            simsize(simsize), boxsize(dx*n),
            dx(dx), x0(x0), y0(y0), z0(z0),
            size(n), size2(n*n), size3(n*n*n)
    {
        // cerr << "Grid ctor " << this <<  endl;
        if(withField) {
            pField = std::make_shared<std::vector<std::complex<T>>>(size3,0);
            pField->shrink_to_fit();
        }
        setKmin();
    }




    Grid(size_t n): simsize(0), boxsize(n),
                    dx(1.0), x0(0),y0(0),z0(0), size(n), size2(n*n), size3(n*n*n)
    {
      // cerr << "Grid ctor size-only" << endl;
      pField=nullptr;
      setKmin();
    }

    virtual ~Grid() {
      // cerr << "~Grid " << this << endl;
    }

    T getWrappedDelta(T x0, T x1) const {
        T result = x0-x1;
        if(result>simsize/2) {
            result-=simsize;
        }
        if(result<-simsize/2) {
            result+=simsize;
        }
        return result;
    }

    static int getRatioAndAssertInteger(T p, T q) {
        const T tolerance = 1e-6;
        T ratio = p/q;
        int rounded_ratio = int(round(ratio));
        assert(abs(T(rounded_ratio)-ratio)<tolerance);
        return rounded_ratio;
    }

    GridPtrType makeProxyGridToMatch(const Grid<T> &target) const
    {
        GridPtrType proxy = std::const_pointer_cast<Grid<T>>(this->shared_from_this());
        if(target.dx>dx) {
            int ratio = getRatioAndAssertInteger(target.dx,dx);
            proxy = std::make_shared<SubSampleGrid<T>>(proxy, ratio);
        } else if(target.dx<dx) {
            int ratio = getRatioAndAssertInteger(dx,target.dx);
            proxy = std::make_shared<SuperSampleGrid<T>>(proxy, ratio);
        }

        if(target.x0!=x0 || target.y0!=y0 || target.z0!=z0 || target.size!=proxy->size) {
            proxy = std::make_shared<SectionOfGrid<T>>(proxy,
                                                       getRatioAndAssertInteger(target.x0-x0, proxy->dx),
                                                       getRatioAndAssertInteger(target.y0-y0, proxy->dx),
                                                       getRatioAndAssertInteger(target.z0-z0, proxy->dx),
                                                       target.size);
        }



        return proxy;
    }

    template<typename TArray>
    void addFieldFromDifferentGrid(const Grid<T> &grid_src, const TArray &pField_x_src, TArray &pField_x_dest) {

        // TODO: make signature above const-correct

        GridPtrType pSourceProxyGrid = grid_src.makeProxyGridToMatch(*this);

        assert(pSourceProxyGrid->fieldIsSuitableSize(pField_x_src));
        assert(fieldIsSuitableSize(pField_x_dest));
        assert(pSourceProxyGrid->size3 == size3);


        #pragma omp parallel for schedule(static)
        for(size_t ind_l=0; ind_l< size3; ind_l++) {
            if (pSourceProxyGrid->isInDomain(ind_l))
                pField_x_dest[ind_l]+=pSourceProxyGrid->getFieldAt(ind_l, pField_x_src);

        }
    }

    void addFieldFromDifferentGrid(const Grid<T> &grid_src) {
        TField & pField_dest = getFieldReal();
        assert(!grid_src.isFieldFourier());
        const TField & pField_src = grid_src.getField();

        addFieldFromDifferentGrid(grid_src, pField_src, pField_dest);

        auto offset_fields_src = grid_src.getOffsetFields();
        auto offset_fields_dest = getOffsetFields();

        if(std::get<0>(offset_fields_src)->size()>0 && std::get<0>(offset_fields_dest)->size()>0) {
            addFieldFromDifferentGrid(grid_src, *std::get<0>(offset_fields_src), *std::get<0>(offset_fields_dest));
            addFieldFromDifferentGrid(grid_src, *std::get<1>(offset_fields_src), *std::get<1>(offset_fields_dest));
            addFieldFromDifferentGrid(grid_src, *std::get<2>(offset_fields_src), *std::get<2>(offset_fields_dest));
        }
    }

    void addFieldFromDifferentGrid(Grid<T> &grid_src) {
        grid_src.getFieldReal();
        addFieldFromDifferentGrid(const_cast<const Grid<T> &>(grid_src));
    }

    virtual void debugInfo(std::ostream& s) const {
        s << "Grid of side " << size << " address " << this;
    }


    virtual void gatherParticleList(std::vector<size_t> & targetArray) const {
      targetArray.insert(targetArray.end(), particleArray.begin(), particleArray.end());
    }

    virtual void distributeParticleList(const std::vector<size_t> & sourceArray) {
      particleArray.insert(particleArray.end(), sourceArray.begin(), sourceArray.end());
    }

    virtual void clearParticleList() {
      particleArray.clear();
    }

    virtual size_t estimateParticleListSize() {
      return particleArray.size();
    }


    virtual bool pointsToGrid(Grid<T> *pOther) {
      return this==pOther;
    }

    bool pointsToAnyGrid(std::vector<std::shared_ptr<Grid<T>>> grids) {
      for(auto g: grids) {
        if(pointsToGrid(g.get()))
          return true;
      }
      return false;
    }

    ///////////////////////////////////////
    //  Field manipulation routines
    ///////////////////////////////////////

    virtual  TField & getFieldFourier()  {
        assert(pField!=nullptr);
        if(!fieldFourier) {
            fft(pField->data(),pField->data(),size,1);
            fieldFourier=true;
        }
        return *pField;
    }

    virtual  TField & getFieldReal()  {
        assert(pField!=nullptr);
        if(fieldFourier) {
            fft(pField->data(),pField->data(),size,-1);
            fieldFourier=false;
        }
        return *pField;
    }

    virtual const TField & getField() const {
        assert(pField!=nullptr);
        return *pField;
    }

    TField & getField() {
        return const_cast<TField &>(const_cast<const Grid *>(this)->getField());
    }

    void applyFilter(const Filter<T> & filter, TField & fieldFourier) {
        assert(fieldFourier.size()==size3);

        #pragma omp parallel for
        for(size_t i=0; i<size3; ++i) {
            fieldFourier[i]*=filter(getAbsK(i));
        }
    }


    virtual bool fieldIsSuitableSize(const TField & field) {
        return field.size()==size3;
    }

    virtual bool fieldIsSuitableSize(const TRealField & field) {
        return field.size()==size3;
    }

    virtual bool isInDomain(size_t i) {
        return i<size3;
    }

    virtual bool isInDomain(int x, int y, int z) {
        return x>=0 && y>=0 && z>=0 && x<size && y<size && z<size;
    }

    virtual complex<T> getFieldAt(size_t i, const TField & field) {
        return field[i];
    }

    virtual T getFieldAt(size_t i, const TRealField & field) {
        return field[i];
    }

    virtual complex<T> getFieldAt(size_t i) {
        return getFieldAt(i, getFieldReal());
    }


    auto getOffsetFields() {
        return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    auto getOffsetFields() const {
        return std::make_tuple(pOff_x, pOff_y, pOff_z);
    }

    virtual bool isFieldFourier()  const {
        return fieldFourier;
    }

    void getParticle(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        getParticleNoWrap(id,x,y,z,vx,vy,vz,cellMassi,eps);
        simWrap(x,y,z);
    }

    virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        getParticleNoOffset(id, x, y, z, vx, vy, vz, cellMassi, eps);
        addCentroidLocation(id,x,y,z);
    }

    virtual T getMass() const {
        return cellMass;
    }

    virtual T getEps() const {
        return dx*0.01075; // <-- arbitrary to coincide with normal UW resolution. TODO: Find a way to make this flexible.
    }

    void getParticleNoOffset(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
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




    virtual std::shared_ptr<Grid<T>> makeScaledMassVersion(T massRatio) {
        return std::make_shared<MassScaledGrid<T>>(this->shared_from_this(), massRatio);
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

    size_t findNextInd(size_t index, const int step[3]) const
    {
        int grid[3];
        std::tie(grid[0],grid[1],grid[2])=getCoordinates(index);

        grid[0]+=step[0];
        grid[1]+=step[1];
        grid[2]+=step[2];

        return this->getIndex(grid); // N.B. does wrapping inside getIndex
    }

    size_t findNextIndNoWrap(size_t index, const int step[3]) const
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
        if(x>(signed) size) x-=size;
        if(y>(signed) size) y-=size;
        if(z>(signed) size) z-=size;
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


    size_t getIndex(int x, int y, int z) const
    {

        size_t size=this->size;

        wrap(x,y,z);

        size_t index=(x*size+y);
        index*=size;
        index+=z;

        return index;

    }

    size_t getIndexNoWrap(int x, int y, int z) const
    {

#ifdef SAFER_SLOWER
        if(x<0 || x>=size || y<0 || y>=size || z<0 || z>=size)
            throw std::runtime_error("Grid index out of range in getIndexNoWrap");
#endif
        return size_t(x*size+y)*size+z;
    }

    size_t getIndex(const int pos[3]) const {
        return getIndex(pos[0], pos[1], pos[2]);
    }

    size_t getIndexNoWrap(const int pos[3]) const {
        return getIndexNoWrap(pos[0], pos[1], pos[2]);
    }

    void getCoordinates(size_t id, int &x, int &y, int &z) const {
        // if(id>=size3) throw std::runtime_error("Index out of range");

        // The following implementation is a little faster than using the
        // modulo operator.
        x = int(id/size2);
        id-=size_t(x)*size2;
        y = int(id/size);
        id-=size_t(y)*size;
        z = int(id);


    }

    tuple<int, int, int> getCoordinates(size_t id) const {
        int x, y, z;

        getCoordinates(id, x,y,z);

        return std::make_tuple(x,y,z);
    }

    void getKCoordinates(size_t id, int &x, int &y, int &z) const {
        getCoordinates(id,x,y,z);
        if(x>(signed) size/2) x=x-size;
        if(y>(signed) size/2) y=y-size;
        if(z>(signed) size/2) z=z-size;
    }

    tuple<int, int, int> getKCoordinates(size_t id) const {
        int x, y, z;

        getKCoordinates(id, x,y,z);

        return std::make_tuple(x,y,z);
    }

    T getKSquared(size_t id) const {
        int x,y,z;
        T res;
        getKCoordinates(id,x,y,z);
        res = x*x+y*y+z*z;
        res*=kMinSquared;
        return res;
    }

    T getAbsK(size_t id) const {
        return sqrt(getKSquared(id));
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

    }

    tuple<T, T, T> getCentroidLocation(size_t id) const {
        T xc,yc,zc;
        getCentroidLocation(id,xc,yc,zc);
        return std::make_tuple(xc,yc,zc);
    }

    size_t getClosestIdNoWrap(T x0c, T y0c, T z0c) {
      int xa=((int) floor((x0c-x0-dx/2)/dx));
      int ya=((int) floor((y0c-y0-dx/2)/dx));
      int za=((int) floor((z0c-z0-dx/2)/dx));
      return getIndexNoWrap(xa,ya,za);

    }

    void appendIdsInCubeToVector(T x0c, T y0c, T z0c, T dxc, vector<size_t> &ids) {
        // return all the grid IDs whose centres lie within the specified cube

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
                    ids.emplace_back(getIndex(x,y,z));
                }
            }
        }
    }

    vector<size_t> getIdsInCube(T x0c, T y0c, T z0c, T dxc) {
        vector<size_t> ids;
        appendIdsInCubeToVector(x0c, y0c, z0c, dxc, ids);
        return ids;
    }




};


template<typename T>
class VirtualGrid : public Grid<T> {
protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;
    GridPtrType pUnderlying;

public:
    VirtualGrid(GridPtrType pUnderlying) :
            Grid<T>(
                    pUnderlying->simsize, pUnderlying->size,
                    pUnderlying->dx, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0, false),
            pUnderlying(pUnderlying) {

    }


    VirtualGrid(GridPtrType pUnderlying, T simsize, T gridsize,
                T dx, T x0, T y0, T z0, bool withField) :
            Grid<T>(simsize, gridsize, dx, x0, y0, z0, withField),
            pUnderlying(pUnderlying) {

    }

    virtual void debugName(std::ostream &s) const {
        s << "VirtualGrid";
    }

    virtual void debugInfo(std::ostream &s) const override {
        debugName(s);
        s << " of side " << this->size << " address " << this << " referencing ";
        pUnderlying->debugInfo(s);
    }

    bool pointsToGrid(Grid<T> *pOther) override {
        return pUnderlying.get() == pOther;
    }

    virtual bool fieldIsSuitableSize(const TRealField &field) override {
        return pUnderlying->fieldIsSuitableSize(field);
    }

    virtual bool fieldIsSuitableSize(const TField &field) override {
        return pUnderlying->fieldIsSuitableSize(field);
    }

    virtual T getFieldAt(size_t i, const TRealField &field) override {
        throw std::runtime_error("getFieldAt is not implemented for this type of VirtualGrid");
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) override {
        throw std::runtime_error("getFieldAt is not implemented for this type of VirtualGrid");
    }

    void gatherParticleList(std::vector<size_t> & targetArray) const override {
      pUnderlying->gatherParticleList(targetArray);
    }

    void distributeParticleList(const std::vector<size_t> & sourceArray) override {
      pUnderlying->distributeParticleList(sourceArray);
    }

    void clearParticleList() override {
      pUnderlying->clearParticleList();
    }

    size_t estimateParticleListSize() override {
      return pUnderlying->estimateParticleListSize();
    }

    virtual TField & getFieldFourier() override {
        return pUnderlying->getFieldFourier();
    }

    virtual TField & getFieldReal() override {
        return pUnderlying->getFieldReal();
    }

    virtual const TField & getField() const override {
        return pUnderlying->getField();
    }

    virtual bool isFieldFourier()  const override  {
        return pUnderlying->isFieldFourier();
    }

    virtual void zeldovich(T hfac, T particlecellMass) override {
        throw std::runtime_error("VirtualGrid - does not contain an actual field in memory");
    }


    virtual T getMass() const override {
        return pUnderlying->getMass();
    }

    void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {
      pUnderlying->getParticleNoWrap(id,x,y,z,vx,vy,vz,cellMassi,eps);
    }

    void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {
      pUnderlying->getParticleFromOffset(x,y,z,vx,vy,vz,cellMassi,eps);
    }

    complex<T> getFieldAt(size_t i) override {
        return this->getFieldAt(i, this->pUnderlying->getFieldReal());
    }

};

template<typename T>
class SuperSampleGrid : public VirtualGrid<T> {
private:
    int factor;
    int factor3;

protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;

public:
    SuperSampleGrid(GridPtrType pUnderlying, int factor):
            VirtualGrid<T>(pUnderlying,
                    pUnderlying->simsize, pUnderlying->size*factor,
                    pUnderlying->dx/factor, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0, false),
            factor(factor)
    {

        factor3=factor*factor*factor;
    }

    virtual T getMass() const override {
        return this->pUnderlying->getMass()/factor3;
    }

    virtual void debugName(std::ostream &s) const override {
        s << "SuperSampleGrid";
    }

    void gatherParticleList(std::vector<size_t> & targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->gatherParticleList(underlyingArray);
      Grid<T>::upscaleParticleList(underlyingArray, targetArray, this->pUnderlying.get(), this);
    }

    void distributeParticleList(const std::vector<size_t> & sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::downscaleParticleList(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->distributeParticleList(targetArray);
    }

    size_t estimateParticleListSize() override {
      return this->pUnderlying->estimateParticleListSize()*factor3;
    }

    virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        this->getCentroidLocation(id,x,y,z);
        getParticleFromOffset(x, y, z, vx, vy, vz, cellMassi, eps);
    }

    void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {
      this->pUnderlying->getParticleFromOffset(x,y,z,vx,vy,vz,cellMassi,eps);
      // adjust mass
      cellMassi/=factor3;
    }

    virtual T getFieldAt(size_t i, const TRealField &field) override {
        T x,y,z;
        this->getCentroidLocation(i,x,y,z);
        return this->pUnderlying->getFieldInterpolated(x,y,z,field);
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) override {
        T x,y,z;
        this->getCentroidLocation(i,x,y,z);
        return this->pUnderlying->getFieldInterpolated(x,y,z,field);
    }

};

template<typename T>
class OffsetGrid : public VirtualGrid<T> {
private:
    T xOffset,yOffset,zOffset;
    int xOffset_i, yOffset_i, zOffset_i;

protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;


public:
    OffsetGrid(GridPtrType pUnderlying, T dx, T dy, T dz):
            VirtualGrid<T>(pUnderlying),
            xOffset(dx), yOffset(dy), zOffset(dz)
    {

    }


    virtual void debugName(std::ostream &s) const override {
        s << "OffsetGrid";
    }

    virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        this->pUnderlying->getParticleNoWrap(id,x,y,z,vx,vy,vz,cellMassi,eps);
        x+=xOffset;
        y+=yOffset;
        z+=zOffset;
        this->simWrap(x,y,z);
    }

    void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {
        x-=xOffset;
        y-=yOffset;
        z-=zOffset;
        this->simWrap(x,y,z);
        this->pUnderlying->getParticleFromOffset(x,y,z,vx,vy,vz,cellMassi,eps);
        x+=xOffset;
        y+=yOffset;
        z+=zOffset;
        this->simWrap(x,y,z);
    }


};


template<typename T>
class SectionOfGrid : public VirtualGrid<T> {
private:
    int xOffset_i, yOffset_i, zOffset_i;
    T xOffset, yOffset, zOffset;

protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;

    size_t mapIndex(size_t sec_id) const {
        int x,y,z;
        this->getCoordinates(sec_id,x,y,z);
        x+=xOffset_i;
        y+=yOffset_i;
        z+=zOffset_i;
        const size_t underlyingSize = this->pUnderlying->size;

        if(x<0 || x>=underlyingSize || y<0 || y>=underlyingSize || z<0 || z>=underlyingSize)
            throw std::out_of_range("Out of range in SectionOfGrid");
        return this->pUnderlying->getIndex(x,y,z);
    }

public:
    SectionOfGrid(GridPtrType pUnderlying, int deltax, int deltay, int deltaz, size_t size):
            VirtualGrid<T>(pUnderlying,
                           pUnderlying->simsize, size,
                           pUnderlying->dx, pUnderlying->x0+deltax*pUnderlying->dx,
                           pUnderlying->y0+deltay*pUnderlying->dx,
                           pUnderlying->z0+deltaz*pUnderlying->dx, false),
            xOffset_i(deltax), yOffset_i(deltay), zOffset_i(deltaz),
            xOffset(deltax*this->dx), yOffset(deltay*this->dx), zOffset(deltaz*this->dx)
    {

    }

    virtual bool isInDomain(size_t i) override {
        int x,y,z;
        this->getCoordinates(i, x,y,z);
        return isInDomain(x,y,z);
    }

    virtual bool isInDomain(int x, int y, int z) override {
        x += xOffset_i;
        y += yOffset_i;
        z += zOffset_i;
        const size_t underlyingSize = this->pUnderlying->size;

        return x >= 0 && y >= 0 && z >= 0 && x < underlyingSize && y < underlyingSize && z < underlyingSize;

    }


    virtual void debugName(std::ostream &s) const override {
        s << "SectionOfGrid";
    }

    virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
        this->pUnderlying->getParticleNoWrap(mapIndex(id),x,y,z,vx,vy,vz,cellMassi,eps);
    }

    void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {

        this->pUnderlying->getParticleFromOffset(x,y,z,vx,vy,vz,cellMassi,eps);

    }

    virtual T getFieldAt(size_t i, const TRealField &field) override {
        return this->pUnderlying->getFieldAt(mapIndex(i), field);
    }

    virtual complex<T> getFieldAt(size_t i, const TField &field) override {
        return this->pUnderlying->getFieldAt(mapIndex(i), field);
    }



};


template<typename T>
class SubSampleGrid : public VirtualGrid<T> {
private:
    int factor;
    int factor3;

protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;

public:
    SubSampleGrid(std::shared_ptr<Grid<T>> pUnderlying, int factor):
            VirtualGrid<T>(pUnderlying,
                    pUnderlying->simsize, pUnderlying->size/factor,
                    pUnderlying->dx*factor, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0, false),
            factor(factor)
    {
        //if(this->pUnderlying->size%factor!=0)
        //  throw std::runtime_error("SubSampleGrid - factor must be a divisor of the original grid size");

        factor3=factor*factor*factor;
    }

    virtual void debugName(std::ostream &s) const override {
        s << "SubSampleGrid";
    }

    void gatherParticleList(std::vector<size_t> & targetArray) const override {
      std::vector<size_t> underlyingArray;
      this->pUnderlying->gatherParticleList(underlyingArray);
      Grid<T>::downscaleParticleList(underlyingArray, targetArray, this->pUnderlying.get(), this);
      // err << "SubSample gatherParticleList - underlying = " << underlyingArray.size() << " transformed = " <<targetArray.size() << endl;
    }

    void distributeParticleList(const std::vector<size_t> & sourceArray) override {
      std::vector<size_t> targetArray;
      Grid<T>::upscaleParticleList(sourceArray, targetArray, this, this->pUnderlying.get());
      this->pUnderlying->distributeParticleList(targetArray);
      // cerr << "SubSample distributeParticleList - source = " << sourceArray.size() << " transformed = " <<targetArray.size() << endl;
    }

    size_t estimateParticleListSize() override {
      return this->pUnderlying->estimateParticleListSize()/factor3;
    }


    virtual T getMass() const override {
        return this->pUnderlying->getMass()*factor3;
    }

    int forEachSubcell(size_t id, std::function<void(size_t)> callback) const {
        int x0,y0,z0;
        this->getCoordinates(id, x0, y0, z0);
        x0*=factor; y0*=factor; z0*=factor;
        auto x1=x0+factor, y1=y0+factor, z1=z0+factor;

        // In case the edge of the last cell in the fine grid (this->pUnderlying)
        // is not aligned with the edge of any cell in the coarse grid (this),
        // we need to be able to do an average over fewer than factor^3 cells.
        //
        // This is a situation we might, in retrospect, wish to avoid. However,
        // since early tests and simulations were conducted without obeying
        // this 'end-alignment' condition, we need to support it.
        if(x1>this->pUnderlying->size) x1=this->pUnderlying->size;
        if(y1>this->pUnderlying->size) y1=this->pUnderlying->size;
        if(z1>this->pUnderlying->size) z1=this->pUnderlying->size;

        int localFactor3 = (x1-x0)*(y1-y0)*(z1-z0);

        for(auto xi=x0; xi<x1; ++xi) {
            for (auto yi=y0; yi<y1; ++yi) {
                for (auto zi=z0; zi<z1; ++zi) {
                    callback(this->pUnderlying->getIndexNoWrap(xi,yi,zi));
                }

            }
        }
        return localFactor3;
    }

    virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
      T xt,yt,zt,vxt,vyt,vzt,cellMassit,epst;

      x=0;
      y=0;
      z=0;
      vx=0;
      vy=0;
      vz=0;
      cellMassi=0;
      eps=0;

      int localFactor3 = forEachSubcell(id, [&](size_t sub_id) {
          this->pUnderlying->getParticleNoWrap(sub_id,
                                               xt,yt,zt,vxt,vyt,vzt,cellMassit,epst);
          x+=xt;
          y+=yt;
          z+=zt;
          vx+=vxt;
          vy+=vyt;
          vz+=vzt;
          eps+=epst;
          cellMassi+=cellMassit;
      });

      // most variables want an average, not a sum:
      x/=localFactor3;
      y/=localFactor3;
      z/=localFactor3;
      vx/=localFactor3;
      vy/=localFactor3;
      vz/=localFactor3;
      eps/=localFactor3;

      // cell mass wants to be a sum over the *entire* cell (even if
      // some subcells are missing):
      cellMassi*=factor3/localFactor3;

    }

    virtual complex<T> getFieldAt(size_t i, const TField & field) {
        complex<T> returnVal(0);
        int localFactor3 = forEachSubcell(i, [this, &returnVal, &field](size_t local_id){
            returnVal+=field[local_id];
        });
        return returnVal/T(localFactor3);
    }

    virtual T getFieldAt(size_t i, const TRealField & field) {
        T returnVal(0);
        int localFactor3 = forEachSubcell(i, [this, &returnVal, &field](size_t local_id){
            returnVal+=field[local_id];
        });
        return returnVal/localFactor3;
    }



    void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {
      this->pUnderlying->getParticleFromOffset(x,y,z,vx,vy,vz,cellMassi,eps);
      cellMassi*=factor3;
    }

};


template<typename T>
class MassScaledGrid : public VirtualGrid<T> {
protected:
    using typename Grid<T>::TField;
    using typename Grid<T>::TRealField;
    using typename Grid<T>::PtrTField;
    using typename Grid<T>::GridPtrType;

    T massFac;

public:
    MassScaledGrid(GridPtrType pUnderlying, T massFac):
            VirtualGrid<T>(pUnderlying,
                    pUnderlying->simsize, pUnderlying->size,
                    pUnderlying->dx, pUnderlying->x0, pUnderlying->y0,
                    pUnderlying->z0, false),
            massFac(massFac)
    {
    }

    virtual void debugName(std::ostream &s) const override {
        s << "MassScaledGrid";
    }

    virtual T getMass() const override {
        return this->pUnderlying->getMass()*massFac;
    }

    virtual void getParticleNoWrap(size_t id, T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const
    {
      this->pUnderlying->getParticleNoWrap(id,x,y,z,vx,vy,vz,cellMassi,eps);
      cellMassi*=massFac;
    }

    void getParticleFromOffset(T &x, T &y, T &z, T &vx, T &vy, T &vz, T &cellMassi, T &eps) const override
    {
      this->pUnderlying->getParticleFromOffset(x,y,z,vx,vy,vz,cellMassi,eps);
      cellMassi*=massFac;
    }

};



#endif

