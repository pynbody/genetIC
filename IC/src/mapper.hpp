//
// mapper.hpp
//
// The idea of the classes in this file are to provide a way to map particle IDs
// onto grid locations. This becomes quite complicated when one has multiple
// grids and potentially things like gas particles and so on. By making
// specific classes to do very specific bits of the mapping, hopefully we keep
// this complexity managable while also having a lot of flexibility to handle
// different set-ups.

#ifndef __MAPPER_HPP
#define __MAPPER_HPP

#include <memory>
#include <typeinfo>
#include "grid.hpp"


// helper function for our debug dumps
void indent(std::ostream& s, int level=0) {
    for(int i=0; i<level; i++) s << "| ";
    s << "+ ";

}

// forward declarations:
template<typename T> class ParticleMapper;
template<typename T> class OneLevelParticleMapper;
template<typename T> class TwoLevelParticleMapper;
template<typename T> class AddGasMapper;

template<typename T>
class MapperIterator {
protected:
    // This is a generic structure used by all subclasses to keep
    // track of where they are in a sequential operation
    size_t i;
    std::vector<std::shared_ptr<MapperIterator<T>>> subIterators;
    std::vector<size_t> extraData;

    std::vector<std::shared_ptr<std::vector<T>>> particleCache;

    const ParticleMapper<T>* pMapper;

    using MapType = ParticleMapper<T>;
    using MapPtrType = std::shared_ptr<ParticleMapper<T>>;
    using GridType = Grid<T>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;
    using DereferenceType = std::pair<ConstGridPtrType, size_t>;

    friend class ParticleMapper<T>;
    friend class OneLevelParticleMapper<T>;
    friend class TwoLevelParticleMapper<T>;
    friend class AddGasMapper<T>;

    MapperIterator(const ParticleMapper<T>* pMapper) : i(0), pMapper(pMapper) {}


public:

    MapperIterator & operator++() {
        pMapper->incrementIterator(this);
        return (*this);
    }

    MapperIterator operator++(int) {
        MapperIterator R = (*this);
        pMapper->incrementIterator(this);
        return R;
    }

    MapperIterator & operator+=(size_t m) {
        pMapper->incrementIteratorBy(this, m);
        return (*this);
    }

    MapperIterator & operator-=(size_t m) {
        pMapper->decrementIteratorBy(this, m);
        return (*this);
    }

    DereferenceType operator*() const {
        return pMapper->dereferenceIterator(this);
    }

    template<typename... Args>
    void getParticle(Args&&... args) {
        const auto q = **this;
        q.first->getParticle(q.second, std::forward<Args>(args)...);
    }

    void precacheParticles() {

    }

    void getNextNParticles(size_t n, std::vector<T> &xAr, std::vector<T> &yAr, std::vector<T> &zAr,
                                     std::vector<T> &vxAr, std::vector<T> &vyAr, std::vector<T> &vzAr) {

    }

    T getMass() {
        const auto q = **this;
        return q.first->getMass();
    }

    std::unique_ptr<DereferenceType> operator->() const {
        return std::unique_ptr<DereferenceType>(new DereferenceType(pMapper->dereferenceIterator(this)));
    }

    friend bool operator==(const MapperIterator<T> & lhs, const MapperIterator<T> & rhs) {
        return (lhs.i == rhs.i) && (lhs.pMapper == rhs.pMapper) ;
    }

    friend bool operator!=(const MapperIterator<T> & lhs, const MapperIterator<T> & rhs) {
        return (lhs.i != rhs.i) || (lhs.pMapper != rhs.pMapper) ;
    }

    void debugInfo(int n=0) const {
        indent(cerr,n);
        cerr << "i=" << i << endl;
        indent(cerr,n);
        cerr << "type=" << typeid(*pMapper).name() << endl;
        for(auto q: extraData) {
            indent(cerr,n);
            cerr << "data=" << q << endl;
        }
        for(auto q: subIterators) {
            indent(cerr,n);
            if(q!=nullptr) {
                cerr << "subiterator: " << endl;
                q->debugInfo(n+1);
            } else {
                cerr << "null subiterator" << endl;
            }
        }
    }


};




template<typename T>
class ParticleMapper {
public:
    using MapType = ParticleMapper<T>;
    using MapPtrType = std::shared_ptr<ParticleMapper<T>>;
    using GridType = Grid<T>;
    using GridPtrType = std::shared_ptr<Grid<T>>;
    using ConstGridPtrType = std::shared_ptr<const Grid<T>>;
    using iterator = MapperIterator<T>;

    friend class MapperIterator<T>;

protected:


    virtual void incrementIterator(iterator *pIterator) const {
        pIterator->i++;
    }

    virtual void incrementIteratorBy(iterator *pIterator,size_t increment) const {
        pIterator->i+=increment;
    }

    virtual void decrementIteratorBy(iterator *pIterator,size_t increment) const {
        throw std::runtime_error("Attempting to reverse in a mapper that does not support random access");
    }

    virtual std::pair<ConstGridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const {
        throw std::runtime_error("There is no grid associated with this particle mapper");
    }

    template<typename... Args>
    void getParticleFromIterator(const iterator *pIterator, Args&&... args) const {

    }

    T getMassFromIterator(const iterator *pIterator) {
        return (*pIterator)->first->getMass();
    }

public:

    virtual void debugInfo(std::ostream& s, int level=0) const {
        indent(s,level);
        s << "Unspecified MapperIterator (shouldn't be here)" << std::endl;
    }

    friend std::ostream& operator<< (std::ostream& stream, const ParticleMapper<T>& I) {
        I.debugInfo(stream);
        return stream;
    }


    virtual size_t size() const {
        return 0;
    }

    virtual size_t size_gas() const {
        return 0;
    }

    virtual size_t size_dm() const {
        return size();
    }

    virtual void interpretParticleList(const std::vector<size_t> & genericParticleArray) {
        throw std::runtime_error("Cannot interpret particles yet; no particle->grid mapper has been set up");
    }

    virtual void interpretParticleList(std::vector<size_t> && genericParticleArray) {
        // Sometimes it's helpful to have move semantics, but in general just call
        // the normal implementation
        this->interpretParticleList(genericParticleArray);
    }

    virtual GridPtrType getCoarsestGrid() {
        throw std::runtime_error("There is no grid associated with this particle mapper");
    }

    virtual iterator begin() const {
        return iterator(this);
    }

    virtual iterator end() const {
        iterator x(this);
        x.i = size();

        return x;
    }

    virtual iterator beginDm() const {
        return begin();
    }

    virtual iterator endDm() const {
        return end();
    }

    virtual iterator beginGas() const {
        return end();
    }

    virtual iterator endGas() const {
        return end();
    }


    virtual bool supportsReverseIterator()  {
        return false;
    }

    virtual MapPtrType addGas(T massRatio, const std::vector<GridPtrType> & toGrids) {
        throw std::runtime_error("Don't know how to add gas in this context");
    }

    virtual MapPtrType superSample(int ratio, const std::vector<GridPtrType> & toGrids) {
        throw std::runtime_error("Don't know how to supersample in this context");
    }





};


template<typename T>
class OneLevelParticleMapper : public ParticleMapper<T> {

private:

    using MapType = ParticleMapper<T>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::ConstGridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;

    GridPtrType pGrid;

protected:

    virtual std::pair<ConstGridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const override {
        return std::make_pair(pGrid,pIterator->i);
    }

    virtual void decrementIteratorBy(iterator *pIterator,size_t increment) const override {
        pIterator->i-=increment;
    }


public:
    virtual void debugInfo(std::ostream& s, int level=0) const override {
        indent(s,level);
        s << "Grid of side " << pGrid->size << std::endl;
    }

    OneLevelParticleMapper(std::shared_ptr<Grid<T>> &pGrid) : pGrid(pGrid) {

    }

    OneLevelParticleMapper(std::shared_ptr<Grid<T>> &&pGrid) : pGrid(pGrid) {

    }

    size_t size() const override {
        return pGrid->size3;
    }

    void interpretParticleList(const std::vector<size_t> & genericParticleArray) override {
        pGrid->particleArray = genericParticleArray;
    }

    void interpretParticleList(std::vector<size_t> && genericParticleArray) override {
        // genericParticleArray is an rvalue - move it instead of copying it
        pGrid->particleArray = std::move(genericParticleArray);
    }

    GridPtrType getCoarsestGrid() override {
        return pGrid;
    }

    bool supportsReverseIterator() override {
        return true;
    }

    MapPtrType addGas(T massRatio, const std::vector<GridPtrType> & toGrids) override
    {
        if(std::find(toGrids.begin(), toGrids.end(), pGrid)!=toGrids.end()) {
            return std::make_shared<OneLevelParticleMapper<T>>(this->pGrid->massSplit(massRatio));
        } else {
            return nullptr;
        }
    }

    MapPtrType superSample(int ratio, const std::vector<GridPtrType> & toGrids) override
    {
        if(std::find(toGrids.begin(), toGrids.end(), pGrid)!=toGrids.end()) {
            auto newGrid = std::make_shared<SuperSampleGrid<T>>(this->pGrid, ratio);
            return std::make_shared<OneLevelParticleMapper<T>>(newGrid);
        } else {
            return std::make_shared<OneLevelParticleMapper<T>>(this->pGrid);
        }
    }





};

template<typename T>
class TwoLevelParticleMapper : public ParticleMapper<T> {
private:

    using MapType = ParticleMapper<T>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;
    using typename MapType::ConstGridPtrType;

    MapPtrType pLevel1;
    MapPtrType pLevel2;

    GridPtrType pGrid1;
    GridPtrType pGrid2;

    unsigned int n_hr_per_lr;

    size_t totalParticles;
    size_t firstLevel2Particle;

    bool skipLevel1;

    std::vector<size_t> mapid(size_t id0) const {
        // Finds all the particles at level2 corresponding to the single
        // particle at level1
        T x0,y0,z0;
        std::tie(x0,y0,z0) = pGrid1->getCentroidLocation(id0);
        return pGrid2->getIdsInCube(x0,y0,z0,pGrid1->dx);
    }


public:
    std::vector<size_t> zoomParticleArray; //< the particles on the coarse grid which we wish to replace with their zooms
protected:
    mutable std::vector<size_t> zoomParticleArrayZoomed; //< the particles on the fine grid that are therefore included

    void syncL2iterator(const size_t & next_zoom, iterator & level2iterator) const {

        if(zoomParticleArrayZoomed[next_zoom]>level2iterator.i) {
            level2iterator+=zoomParticleArrayZoomed[next_zoom]-level2iterator.i;
        } else {
            level2iterator-=level2iterator.i-zoomParticleArrayZoomed[next_zoom];
        }

        assert(level2iterator.i==zoomParticleArrayZoomed[next_zoom]);
    }

    virtual void incrementIterator(iterator *pIterator) const override {

        // set up some helpful shortcut references
        auto & extraData = pIterator->extraData;
        iterator & level1iterator = *(pIterator->subIterators[0]);
        iterator & level2iterator = *(pIterator->subIterators[1]);
        size_t & i = pIterator->i;

        size_t & next_zoom = extraData[0];
        size_t & next_zoom_index = extraData[1];

        // increment the boring-old-counter!
        i++;

        // now work out what it actually points to...

        if(i>=firstLevel2Particle) {
            // do zoom particles

            if(i==firstLevel2Particle)
                next_zoom=0;
            else
                next_zoom++;

            syncL2iterator(next_zoom,level2iterator);


        } else {
            // do normal particles
            //
            // N.B. should never reach here if level1iterator is NULL

            assert((&level1iterator)!=nullptr);

            ++level1iterator;

            while(level1iterator.i==next_zoom_index) {
                // we DON'T want to return this particle.... it will be 'zoomed'
                // later on... ignore it
                ++level1iterator;

                // load the next zoom particle
                ++next_zoom;
                if(next_zoom<zoomParticleArray.size())
                    next_zoom_index = zoomParticleArray[next_zoom];
                else
                    next_zoom_index = size()+1; // i.e. there isn't a next zoom index!

                // N.B. now it's possible we our 'next' particle is also to be zoomed on,
                // so the while loop goes back round to account for this
            }



        }

    }


    virtual void incrementIteratorBy(iterator *pIterator,size_t increment) const {
        // could be optimized:
        for(size_t i =0; i<increment; i++)
            incrementIterator(pIterator);

    }

    virtual std::pair<ConstGridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const override {
        if(pIterator->i>=firstLevel2Particle)
            return **(pIterator->subIterators[1]);
        else
            return **(pIterator->subIterators[0]);
    }

    void calculateHiresParticleList() const {

        for(size_t i=0; i<zoomParticleArray.size(); i++) {
            // get the list of zoomed particles corresponding to this one
            std::vector<size_t> hr_particles = mapid(zoomParticleArray[i]);
            assert(hr_particles.size()==n_hr_per_lr);
            zoomParticleArrayZoomed.insert(zoomParticleArrayZoomed.end(),
                                           hr_particles.begin(), hr_particles.end());

        }

        if(!pLevel2->supportsReverseIterator()) {
            // underlying map can't cope with particles being out of order - sort them
            std::sort(zoomParticleArrayZoomed.begin(), zoomParticleArrayZoomed.end() );
        }

    }





public:

    virtual iterator begin() const override {
        iterator x(this);
        x.extraData.push_back(0); // current position in zoom ID list
        x.extraData.push_back(zoomParticleArray[0]); // next entry in zoom ID list

        if(!skipLevel1)
            x.subIterators.emplace_back(new iterator(pLevel1->begin()));
        else
            x.subIterators.emplace_back(nullptr);

        x.subIterators.emplace_back(new iterator(pLevel2->begin()));

        syncL2iterator(0,*(x.subIterators.back()));

        return x;
    }


    TwoLevelParticleMapper(  MapPtrType & pLevel1,
                             MapPtrType & pLevel2,
                             const std::vector<size_t> & zoomParticles,
                             int n_hr_per_lr,
                             bool skipLevel1=false) :
             pLevel1(pLevel1), pLevel2(pLevel2),
             pGrid1(pLevel1->getCoarsestGrid()),
             pGrid2(pLevel2->getCoarsestGrid()),
             n_hr_per_lr(n_hr_per_lr),
             skipLevel1(skipLevel1),
             zoomParticleArray(zoomParticles)
    {

        std::sort(zoomParticleArray.begin(), zoomParticleArray.end() );

        totalParticles = pLevel1->size()+(n_hr_per_lr-1)*zoomParticleArray.size();

        firstLevel2Particle = pLevel1->size()-zoomParticleArray.size();

        if(skipLevel1) {
            firstLevel2Particle = 0;
            totalParticles = n_hr_per_lr*zoomParticleArray.size();
        }

        assert(pLevel1->size_gas()==0);
        assert(pLevel2->size_gas()==0);

        calculateHiresParticleList();


    }

    virtual void debugInfo(std::ostream& s, int level=0) const override {
        indent(s,level);
        s << "TwoLevelParticleMapper, n_hr_per_lr=" << n_hr_per_lr << ", firstLevel2Particle=" << firstLevel2Particle << std::endl;
        indent(s,level);
        s << "                      , zoom.size=" << zoomParticleArray.size() << ", zoomed.size=" << zoomParticleArrayZoomed.size() << std::endl;
        if(skipLevel1) {
            indent(s,level);
            s << "low-res part will be skipped but notionally is:" << endl;
        }
        pLevel1->debugInfo(s,level+1);
        indent(s,level);
        s << "TwoLevelParticleMapper continues with high-res particles:" << std::endl;
        pLevel2->debugInfo(s,level+1);
        indent(s,level);
        s << "TwoLevelParticleMapper ends" << std::endl;
    }

    void interpretParticleList(const std::vector<size_t> & genericParticleArray) override {

        std::vector<size_t> grid1particles;
        std::vector<size_t> grid2particles;

        if(pLevel1->size()<zoomParticleArray.size())
            throw std::runtime_error("Zoom particle list is longer than the grid it refers to");

        size_t i0 = pLevel1->size()-zoomParticleArray.size();
        size_t lr_particle;


        for(auto i=genericParticleArray.begin(); i!=genericParticleArray.end(); ++i) {
            if((*i)<i0)
                throw std::runtime_error("Constraining particle is in low-res region - not permitted");

            if( ((*i)-i0)/n_hr_per_lr >= zoomParticleArray.size() )
                throw std::runtime_error("Particle ID out of range");

            // find the low-res particle. Note that we push it onto the list
            // even if it's already on there to get the weighting right and
            // prevent the need for a costly search
            lr_particle = zoomParticleArray[((*i)-i0)/n_hr_per_lr];

            grid1particles.push_back(lr_particle);

            // get all the HR particles
            std::vector<size_t> hr_particles = mapid(lr_particle);
            assert(hr_particles.size()==n_hr_per_lr);

            // work out which of these this particle must be and push it onto
            // the HR list
            int offset = ((*i)-i0)%n_hr_per_lr;
            grid2particles.push_back(hr_particles[offset]);

        }

        // Get the grids to interpret their own particles. This almost
        // certainly just involves making a copy of our list. In fact, we're
        // done with our list so they might as well steal the data.
        pLevel1->interpretParticleList(std::move(grid1particles));
        pLevel2->interpretParticleList(std::move(grid2particles));

    }

    virtual GridPtrType getCoarsestGrid() override {
        return pGrid1;
    }

    virtual size_t size() const override {
        return totalParticles;
    }

    MapPtrType addGas(T massRatio, const std::vector<GridPtrType> & toGrids) override
    {
        bool newskip = skipLevel1;

        auto gsub1 = pLevel1->addGas(massRatio, toGrids);
        auto gsub2 = pLevel2->addGas(massRatio, toGrids);

        if(gsub1==nullptr) {
            gsub1 = pLevel1;
            newskip = true;
        }

        if(gsub2!=nullptr)
            return std::make_shared<TwoLevelParticleMapper<T>>(
                gsub1, gsub2, zoomParticleArray,
                n_hr_per_lr, newskip);
        else
            return nullptr;
    }

    MapPtrType superSample(int ratio, const std::vector<GridPtrType> & toGrids) override
    {

        auto ssub1 = pLevel1->superSample(ratio, toGrids);
        int upgrade1 = ssub1->size()/pLevel1->size();

        auto ssub2 = pLevel2->superSample(ratio, toGrids);
        int upgrade2 = ssub2->size()/pLevel2->size();

        return std::make_shared<TwoLevelParticleMapper<T>>(
            ssub1, ssub2, zoomParticleArray,
            n_hr_per_lr*upgrade2/upgrade1,
            skipLevel1);
    }

};

/*
template<typename T>
class SubMapper : public ParticleMapper<T>
{
public:
    using MapType = ParticleMapper<T>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;
    using typename MapType::ConstGridPtrType;

protected:
    size_t start;
    size_t finish;
    MapPtrType pUnderlying;

    virtual void incrementIterator(iterator *pIterator) const {
        pIterator->i++;
        (*(pIterator->subIterators[0]))++;
    }


    virtual std::pair<ConstGridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const {
        return **(pIterator->subIterators[0]);
    }


public:

    virtual size_t size() const {
        return finish-start;
    }


    SubMapper( MapPtrType & pUnderlying, size_t start, size_t finish) :
            pUnderlying(pUnderlying), start(start), finish(finish)
    {
        assert(pUnderlying->size_gas()==0);
    }


    virtual void interpretParticleList(const std::vector<size_t> & genericParticleArray) {
        std::vector<size_t> result;

        for(auto i: genericParticleArray) {
            if(i+start>finish)
            throw std::runtime_error("Out of range!");
            result.push_back(i+start);
        }
        pUnderlying->interpretParticleList(result);

    }

    virtual iterator begin() const {
        iterator i(this);
        i.subIterators.emplace_back(new iterator(pUnderlying->begin()));
        *(i.subIterators[0])+=start;
        return i;
    }

};
*/



template<typename T>
class AddGasMapper : public ParticleMapper<T>
{
public:
    using MapType = ParticleMapper<T>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;
    using typename MapType::ConstGridPtrType;

protected:
    MapPtrType firstMap;
    MapPtrType secondMap;
    bool gasFirst;
    size_t dm0;
    size_t nFirst, nSecond;


    virtual void incrementIterator(iterator *pIterator) const {

        if(pIterator->i>=nFirst)
            (*(pIterator->subIterators[1]))++;
        else
            (*(pIterator->subIterators[0]))++;
        pIterator->i++;

    }


    virtual std::pair<ConstGridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const {
        if(pIterator->i>=nFirst)
            return **(pIterator->subIterators[1]);
        else
            return **(pIterator->subIterators[0]);
    }


public:

    virtual void debugInfo(std::ostream& s, int level=0) const override {
        indent(s,level);
        s << "AddGasMapper";
        if(gasFirst)
            s << ", gas first:" << std::endl;
        else
            s << ", DM first:" << std::endl;

        firstMap->debugInfo(s,level+1);

        indent(s,level);

        s << "AddGasMapper";
        if(gasFirst)
            s << ", DM second:" << std::endl;
        else
            s << ", gas second:" << std::endl;

        secondMap->debugInfo(s,level+1);

        indent(s,level);
        s << "AddGasMapper ends" << std::endl;
    }

    virtual size_t size() const {
        return nFirst+nSecond;
    }

    virtual size_t size_gas() const {
        return gasFirst?nFirst:nSecond;
    }

    virtual size_t size_dm() const {
        return gasFirst?nSecond:nFirst;
    }

    AddGasMapper( MapPtrType & pFirst, MapPtrType &pSecond, bool gasFirst=true ) :
    firstMap(pFirst), secondMap(pSecond), gasFirst(gasFirst), nFirst(pFirst->size()), nSecond(pSecond->size())
    {
        assert(pFirst->size_gas()==0);
        assert(pSecond->size_gas()==0);
    };


    virtual void interpretParticleList(const std::vector<size_t> & genericParticleArray) {
        std::vector<size_t> first;
        std::vector<size_t> second;

        for(auto i: genericParticleArray) {
            if(i<nFirst)
                first.push_back(i);
            else
                second.push_back(i-nFirst);
        }

        firstMap->interpretParticleList(first);
        secondMap->interpretParticleList(second);

    }

    virtual iterator begin() const override {
        iterator i(this);
        i.subIterators.emplace_back(new iterator(firstMap->begin()));
        i.subIterators.emplace_back(new iterator(secondMap->begin()));
        return i;
    }

    virtual iterator beginDm() const override {
        return (gasFirst?secondMap:firstMap)->begin();
    }

    virtual iterator endDm() const override {
        return (gasFirst?secondMap:firstMap)->end();
    }

    virtual iterator beginGas() const override {
        return (gasFirst?firstMap:secondMap)->begin();
    }

    virtual iterator endGas() const override {
        return (gasFirst?firstMap:secondMap)->end();
    }

    MapPtrType superSample(int ratio, const std::vector<GridPtrType> & toGrids) override
    {
        auto ssub1 = firstMap->superSample(ratio, toGrids);
        auto ssub2 = secondMap->superSample(ratio, toGrids);
        cerr << "GAS CONSTRUCT SUPERSAMPLE " << ssub1->size() <<  " " << ssub2->size() << endl;
        return std::make_shared<AddGasMapper<T>>(
            ssub1, ssub2, gasFirst);
    }


};

template class MapperIterator<double>;

#endif // __MAPPER_HPP
