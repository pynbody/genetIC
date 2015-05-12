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
#include "grid.hpp"


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

    MapperIterator(const ParticleMapper<T>* pMapper) : pMapper(pMapper), i(0) {}

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

    MapperIterator & operator+=(size_t i) {
        pMapper->incrementIteratorBy(this, i);
        return (*this);
    }

    MapperIterator & operator-=(size_t i) {
        pMapper->decrementIteratorBy(this, i);
        return (*this);
    }

    DereferenceType operator*() const {
        return pMapper->dereferenceIterator(this);
    }

    template<typename... Args>
    void getParticle(Args&&... args) {
        pMapper->getParticleFromIterator(this, std::forward<Args>(args)...);
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
        const auto q = **pIterator;
        q.first->get_particle(q.second, std::forward<Args>(args)...);
    }



public:

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

    int n_hr_per_lr;

    size_t totalParticles;
    size_t firstLevel2Particle;

    friend class TwoLevelParticleMapper<T>;

    std::vector<size_t> mapid(size_t id0) const {
        // Finds all the particles at level2 corresponding to the single
        // particle at level1
        T x0,y0,z0;
        std::tie(x0,y0,z0) = pGrid1->get_centroid_location(id0);
        return pGrid2->get_ids_in_cube(x0,y0,z0,pGrid1->dx);
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

        // the last 'low-res' particle is just before this address:
        size_t i_write = pLevel1->size();-zoomParticleArray.size();


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

        if(pLevel1!=nullptr)
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
                             const std::vector<size_t> & zoomParticlesZoomed,
                             int n_hr_per_lr) :
                             pLevel1(pLevel1), pLevel2(pLevel2),
                             zoomParticleArray(zoomParticles),
                             n_hr_per_lr(n_hr_per_lr),
                             zoomParticleArrayZoomed(zoomParticlesZoomed)
    {

        // Protected constructor that can deal with the case where we ONLY want to
        // return the zoomed particles from level2

        if(pLevel2==nullptr)
            throw std::runtime_error("Cannot create two level particle mapper with nothing on the second level");

        pGrid2 = pLevel2->getCoarsestGrid();

        if(pLevel1!=nullptr) {
            totalParticles = pLevel1->size()+(n_hr_per_lr-1)*zoomParticleArray.size();
            firstLevel2Particle = pLevel1->size()-zoomParticleArray.size();
            pGrid1 = pLevel1->getCoarsestGrid();
            assert(pLevel1->size_gas()==0);
        } else {
            totalParticles = zoomParticlesZoomed.size();
            firstLevel2Particle=0;
        }

        assert(pLevel2->size_gas()==0);


    }


    TwoLevelParticleMapper(  MapPtrType & pLevel1,
                             MapPtrType & pLevel2,
                             const std::vector<size_t> & zoomParticles,
                             int n_hr_per_lr) :
                             pLevel1(pLevel1), pLevel2(pLevel2),
                             zoomParticleArray(zoomParticles),
                             n_hr_per_lr(n_hr_per_lr),
                             pGrid1(pLevel1->getCoarsestGrid()),
                             pGrid2(pLevel2->getCoarsestGrid())
    {

        std::sort(zoomParticleArray.begin(), zoomParticleArray.end() );

        totalParticles = pLevel1->size()+(n_hr_per_lr-1)*zoomParticleArray.size();

        firstLevel2Particle = pLevel1->size()-zoomParticleArray.size();

        size_t level1_size = pLevel1->size();

        assert(pLevel1->size_gas()==0);
        assert(pLevel2->size_gas()==0);

        calculateHiresParticleList();

    }

    void interpretParticleList(const std::vector<size_t> & genericParticleArray) override {

        std::vector<size_t> grid1particles;
        std::vector<size_t> grid2particles;

        long i0 = pLevel1->size()-zoomParticleArray.size();
        long lr_particle;


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
        auto gsub1 = pLevel1->addGas(massRatio, toGrids);
        auto gsub2 = pLevel2->addGas(massRatio, toGrids);
        return std::make_shared<TwoLevelParticleMapper<T>>(
            gsub1, gsub2, zoomParticleArray, zoomParticleArrayZoomed,
            n_hr_per_lr);
    }

    MapPtrType superSample(int ratio, const std::vector<GridPtrType> & toGrids) override
    {
        auto ssub1 = pLevel1->superSample(ratio, toGrids);
        auto ssub2 = pLevel2->superSample(ratio, toGrids);
        return std::make_shared<TwoLevelParticleMapper<T>>(
            ssub1, ssub2, zoomParticleArray,
            n_hr_per_lr*(ssub2->size()/pLevel2->size())/(ssub1->size()/pLevel1->size())
        );
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
    size_t nSecond;
    size_t nFirst;

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
    firstMap(pFirst), secondMap(pSecond), nFirst(pFirst->size()), nSecond(pSecond->size()), gasFirst(gasFirst)
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
        return std::make_shared<AddGasMapper<T>>(
            firstMap, secondMap, gasFirst);
    }


};

#endif // __MAPPER_HPP
