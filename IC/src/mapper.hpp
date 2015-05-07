//
// mapper.hpp
//
// The idea of the classes in this file are to provide a way to map particle IDs
// onto grid locations. This becomes quite complicated when one has multiple
// grids and potentially things like gas particles and so on. By making
// specific classes to do very specific bits of the mapping, hopefully we keep
// this complexity managable while also having a lot of flexibility to handle
// different set-ups.


#include <memory>


// forward declarations:
template<typename T> class ParticleMapper;
template<typename T> class OneLevelParticleMapper;
template<typename T> class TwoLevelParticleMapper;

template<typename T>
class MapperIterator {
protected:
    // This is a generic structure used by all subclasses to keep
    // track of where they are in a sequential operation
    size_t i;
    std::vector<std::shared_ptr<MapperIterator<T>>> subIterators;
    std::vector<size_t> extraData;
    ParticleMapper<T>* pMapper;

    using MapType = ParticleMapper<T>;
    using MapPtrType = std::shared_ptr<ParticleMapper<T>>;
    using GridType = Grid<T>;
    using GridPtrType = std::shared_ptr<Grid<T>>;

    friend class ParticleMapper<T>;
    friend class OneLevelParticleMapper<T>;
    friend class TwoLevelParticleMapper<T>;

    MapperIterator(ParticleMapper<T>* pMapper) : pMapper(pMapper), i(0) {}

public:

    const MapperIterator & operator++() {
        pMapper->incrementIterator(this);
        return (*this);
    }

    const MapperIterator & operator+=(size_t i) {
        pMapper->incrementIteratorBy(this, i);
        return (*this);
    }

    std::pair<GridPtrType, size_t> operator*() const {
        return pMapper->dereferenceIterator(this);
    }

    friend bool operator==(const MapperIterator<T> & lhs, const MapperIterator<T> & rhs) {
        return (lhs.i == rhs.i) && (lhs.pMapper == rhs.pMapper) ;
    }

    friend bool operator!=(const MapperIterator<T> & lhs, const MapperIterator<T> & rhs) {
        return (lhs.i != rhs.i) || (lhs.pMapper != rhs.pMapper) ;
    }


};

template<typename MyFloat>
class ParticleMapper {
public:
    using MapType = ParticleMapper<MyFloat>;
    using MapPtrType = std::shared_ptr<ParticleMapper<MyFloat>>;
    using GridType = Grid<MyFloat>;
    using GridPtrType = std::shared_ptr<Grid<MyFloat>>;
    using iterator = MapperIterator<MyFloat>;

protected:


    virtual void incrementIterator(iterator *pIterator) const {
        pIterator->i++;
    }

    virtual void incrementIteratorBy(iterator *pIterator,size_t increment) const {
        pIterator->i+=increment;
    }

    virtual std::pair<GridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const {
        throw std::runtime_error("There is no grid associated with this particle mapper");
    }

public:

    virtual size_t size() {
        return 0;
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

    virtual void gatherParticlesOnCoarsestGrid() {
        throw std::runtime_error("There is no grid associated with this particle mapper");
    }

    virtual iterator begin() {
        iterator x(this);
        return x;
    }

    virtual iterator end() {
        iterator x(this);
        x.i = size();
        return x;
    }




};


template<typename MyFloat>
class OneLevelParticleMapper : public ParticleMapper<MyFloat> {

private:

    using MapType = ParticleMapper<MyFloat>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;

    GridPtrType pGrid;

protected:

    virtual std::pair<GridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const override {
        return std::make_pair(pGrid,pIterator->i);
    }

public:
    OneLevelParticleMapper(std::shared_ptr<Grid<MyFloat>> &pGrid) : pGrid(pGrid) {

    }

    OneLevelParticleMapper(std::shared_ptr<Grid<MyFloat>> &&pGrid) : pGrid(pGrid) {

    }

    size_t size() override {
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

    void gatherParticlesOnCoarsestGrid() override {
        // null operation
    }





};

template<typename MyFloat>
class TwoLevelParticleMapper : public ParticleMapper<MyFloat> {
private:

    using MapType = ParticleMapper<MyFloat>;
    using typename MapType::MapPtrType;
    using typename MapType::GridPtrType;
    using typename MapType::GridType;
    using typename MapType::iterator;


    MapPtrType pLevel1;
    MapPtrType pLevel2;

    GridPtrType pGrid1;
    GridPtrType pGrid2;

    int n_hr_per_lr;

    size_t totalParticles;
    size_t firstLevel2Particle;

    std::vector<size_t> mapid(size_t id0) {
        // Finds all the particles at level2 corresponding to the single
        // particle at level1
        MyFloat x0,y0,z0;
        std::tie(x0,y0,z0) = pGrid1->get_centroid_location(id0);
        return pGrid2->get_ids_in_cube(x0,y0,z0,pGrid1->dx);
    }

    void deleteParticles(std::vector<MyFloat*> A) {

        // delete the particles in the 'zoom' region.
        // This involves just moving everything backwards
        // in steps.

        size_t write_ptcl=0;
        size_t next_zoom=zoomParticleArray[0];
        size_t next_zoom_i=0;

        size_t level1_size = pLevel1->size();

        for(size_t read_ptcl=0; read_ptcl<level1_size; read_ptcl++) {
            if(read_ptcl!=next_zoom) {
                if(read_ptcl!=write_ptcl) {
                    for(auto ar=A.begin(); ar!=A.end(); ++ar) {
                        (*ar)[write_ptcl]=(*ar)[read_ptcl];
                    }
                }
                write_ptcl++;
            } else {
                // we've just hit a 'zoom' particle. Skip over it.
                // Keep track of the next zoom particle
                if(next_zoom_i+1<zoomParticleArray.size()) {
                    next_zoom_i++;
                    next_zoom = zoomParticleArray[next_zoom_i];
                } else {
                    // no more zoom particles
                    next_zoom =-1;
                }
            }
        }
    }

    void insertParticles(std::vector<MyFloat*> A, std::vector<MyFloat*> B) {
        // insert the particles from the zoom region

        size_t level1_size = pLevel1->size();

        // the last 'low-res' particle is just before this address:
        size_t i_write = level1_size-zoomParticleArray.size();


        for(size_t i=0; i<zoomParticleArray.size(); i++) {
            // get the list of zoomed particles corresponding to this one
            std::vector<size_t> hr_particles = mapid(zoomParticleArray[i]);

            assert(hr_particles.size()==n_hr_per_lr);

            for(auto i=hr_particles.begin(); i!=hr_particles.end(); i++) {
                for(auto ar_to=A.begin(), ar_from=B.begin();
                    ar_to!=A.end() && ar_from!=B.end();)
                {
                    (*ar_to)[i_write] = (*ar_from)[(*i)];
                    ++ar_from;
                    ++ar_to;
                }
                ++i_write;
            }
        }
        assert(i_write==size());
    }

protected:
    std::vector<size_t> zoomParticleArray; //< the particles on the coarse grid which we wish to replace with their zooms
    std::vector<size_t> zoomParticleArrayZoomed; //< the particles on the fine grid that are therefore included

    virtual void incrementIterator(iterator *pIterator) const override {

        // set up some helpful shortcut references
        auto & extraData = pIterator->extraData;
        iterator & level1iterator = pIterator->subIterators[0];
        iterator & level2iterator = pIterator->subIterators[1];
        size_t & i = pIterator->i;

        if(i>=firstLevel2Particle) {
            // do zoom particles
            if(i==firstLevel2Particle) {
                // initialize the internal iterator structure for the zoom particles
                extraData.clear();
                extraData.push_back(0); // three zeros corresponding to quantities defined below
                extraData.push_back(0);
                extraData.push_back(n_hr_per_lr);
            }

            size_t & zoom_index = extraData[0];
            size_t & buffer_index = extraData[1];
            size_t & buffer_length = extraData[2];

            if(buffer_index==buffer_length) {
                // we have run out of particles to pump out.. get another one
                pIterator->particleArray = mapid(zoomParticleArray[zoom_index]);
                zoom_index++;
                buffer_index = 0;
                assert(zoomParticleArray.size()==n_hr_per_lr);
            }

            level2iterator+=pIterator->particleArray[buffer_index]-level2iterator.i;


        } else {
            // do normal particles
            //
            // here is what the extra data means:
            size_t & next_zoom = extraData[0];
            size_t & next_zoom_index = extraData[1];

            while(level1iterator->i==next_zoom_index) {
                // we DON'T want to return this particle.... it will be 'zoomed'
                // later on... ignore it
                level1iterator++;

                // load the next zoom particle
                next_zoom++;
                if(next_zoom<zoomParticleArray.size())
                    next_zoom_index = zoomParticleArray[next_zoom];
                else
                    next_zoom_index = size()+1; // i.e. there isn't a next zoom index!

                // N.B. now it's possible we our 'next' particle is also to be zoomed on,
                // so the while loop goes back round to account for this
            }



        }


        // increment the boring-old-counter!
        i++;



    }


    virtual void incrementIteratorBy(iterator *pIterator,size_t increment) const {
        // could be optimized:
        for(size_t i =0; i<increment; i++)
            incrementIterator(pIterator);

    }

    virtual std::pair<GridPtrType, size_t> dereferenceIterator(const iterator *pIterator) const override {

    }

public:

    virtual iterator begin() override {
        iterator x(this);
        x.extraData.push_back(0); // current position in zoom ID list
        return x;
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

        // the last 'low-res' particle is just before this address:
        size_t i_write = level1_size-zoomParticleArray.size();


        for(size_t i=0; i<zoomParticleArray.size(); i++) {
            // get the list of zoomed particles corresponding to this one
            std::vector<size_t> hr_particles = mapid(zoomParticleArray[i]);
            assert(hr_particles.size()==n_hr_per_lr);
            zoomParticleArrayZoomed.insert(zoomParticleArrayZoomed.end(),
                                           hr_particles.begin(), hr_particles.end());

        }

        std::sort(zoomParticleArrayZoomed.begin(), zoomParticleArrayZoomed.end() );


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

    virtual size_t size() override {
        return totalParticles;
    }


    void gatherParticlesOnCoarsestGrid() override {
        pLevel1->gatherParticlesOnCoarsestGrid();
        pLevel2->gatherParticlesOnCoarsestGrid();
        deleteParticles(pGrid1->particleProperties);
        insertParticles(pGrid1->particleProperties,pGrid2->particleProperties);
    }



};
