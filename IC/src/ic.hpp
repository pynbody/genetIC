#ifndef _IC_HPP_INCLUDED
#define _IC_HPP_INCLUDED

#include <string>
#include <tuple>
#include <cassert>
#include <functional>
#include <algorithm>
#include <memory>
#include <limits>
#include <iostream>
#include <list>

#include "numpy.hpp"
#include "onelevelconstraint.hpp"
#include "constraintapplicator.hpp"
#include "filesystem.h"
#include "randomfieldgenerator.hpp"
#include "MultiLevelConstraintGenerator.h"

#define for_each_level(level) for(int level=0; level<2 && n[level]>0; level++)

using namespace std;

template<typename MyFloat>
class DummyIC;




template<typename MyFloat>
class IC {
protected:

    using GridPtrType = std::shared_ptr<Grid<MyFloat>>;
    friend class DummyIC<MyFloat>;

    CosmologicalParameters<MyFloat> cosmology;
    MultiLevelFieldManager<MyFloat> fieldManager;
    ConstraintApplicator<MyFloat> constraintApplicator;
    MultiLevelConstraintGenerator<MyFloat> constraintGenerator;
    RandomFieldGenerator<MyFloat> randomFieldGenerator;
    CAMB<MyFloat> spectrum;

    // Everything about the grids:
    MyFloat boxlen[2], dx[2];      // the box length of each grid
    int n[2];                      // number of grid divisions along one axis

    int supersample, subsample;               // DM supersampling to perform on zoom grid, and subsampling on base grid

    std::vector<GridPtrType> pGrid;       // the objects that help us relate points on the grid


    MyFloat x_off[2], y_off[2], z_off[2]; // x,y,z offsets for subgrids

    MyFloat xOffOutput, yOffOutput, zOffOutput;

    int zoomfac; // the zoom factor, i.e. the ratio dx/dx[1]



    int out, gadgetformat;

    size_t nPartLevel[2];

    string incamb, indir, inname, base;


    bool prepared;




    std::vector<size_t> genericParticleArray;
    std::vector<size_t> zoomParticleArray;


    MyFloat x0, y0, z0;

    shared_ptr<ParticleMapper<MyFloat>> pMapper;
    shared_ptr<ParticleMapper<MyFloat>> pInputMapper;

    using RefFieldType = decltype(pGrid[0]->getField());
    using FieldType = std::remove_reference_t<decltype(pGrid[0]->getField())>;

    ClassDispatch<IC<MyFloat>, void> &interpreter;


public:
    IC(ClassDispatch<IC<MyFloat>, void> &interpreter) : interpreter(interpreter),
                                                        pMapper(new ParticleMapper<MyFloat>()),
                                                        constraintApplicator(&fieldManager),
                                                        randomFieldGenerator(fieldManager),
                                                        constraintGenerator(fieldManager, cosmology)
    {
        pInputMapper = nullptr;
        cosmology.hubble =0.701;   // old default
        cosmology.OmegaBaryons0 =-1.0;
        cosmology.ns = 0.96;      // old default
        n[0]=-1;
        n[1]=-1; // no subgrid by default
        boxlen[0]=-1;
        x_off[0]=y_off[0]=z_off[0]=0;
        x_off[1]=y_off[1]=z_off[1]=0;
        prepared = false;
        supersample = 1;
        subsample = 1;
        xOffOutput = 0;
        yOffOutput = 0;
        zOffOutput = 0;
    }

    ~IC() {

    }

    void setOmegaM0(MyFloat in) {
        cosmology.OmegaM0 =in;
    }

    void setOmegaB0(MyFloat in) {
        cosmology.OmegaBaryons0 =in;
        // now that we have gas, mapper may have changed:
        updateParticleMapper();
    }

    void setOmegaLambda0(MyFloat in) {
        cosmology.OmegaLambda0 =in;
    }

    void setHubble(MyFloat in) {
        cosmology.hubble =in;
    }

    void offsetOutput(MyFloat x, MyFloat y, MyFloat z) {
        xOffOutput = x;
        yOffOutput = y;
        zOffOutput = z;
        updateParticleMapper();
    }

    void setSigma8(MyFloat in) {
        cosmology.sigma8 = in;
    }

    void setSupersample(int in) {
        supersample = in;
        updateParticleMapper();
    }

    void setSubsample(int in) {
        subsample = in;
        updateParticleMapper();
    }

    void setBoxLen(MyFloat in) {
        boxlen[0] = in;
        initGrid();
    }

    void setZ0(MyFloat in) {
        cosmology.redshift = in;
        cosmology.scalefactor =1./(cosmology.redshift +1.);
    }

    void setn(int in) {
        n[0] = in;
        initGrid();
    }

    virtual void initGrid(unsigned int level=0) {

        if(n[level]<0 || boxlen[level]<0)
            return;

        nPartLevel[level] = ((size_t) n[level]*n[level])*n[level];
        dx[level] = boxlen[level]/n[level];

        if(pGrid.size()!=level)
            throw std::runtime_error("Trying to re-initialize a grid level");

        pGrid.push_back(std::make_shared<Grid<MyFloat>>(boxlen[0], n[level], dx[level],
                                             x_off[level], y_off[level], z_off[level]));


        updateParticleMapper();
        updateFieldManager();

    }

    void setns(MyFloat in) {
        cosmology.ns = in;
    }

    void setn2(int in) {
        n[1] = in;
        nPartLevel[1] = ((size_t)n[1]*n[1])*n[1];
    }

    void setZoom(int in) {
        // open a subgrid which is the specified factor smaller than
        // the parent grid
        boxlen[1] = boxlen[0]/in;
        zoomfac = in;
    }

    void setZoomParticles(string fname) {
      appendParticleIdFile(fname);
      doZoom();
    }

    void doZoom() {
        if(n[1]==0)
            throw(std::runtime_error("Set n2 before specifying the zoom particles"));

        if(boxlen[1]==0)
            throw(std::runtime_error("Set the zoom factor before specifying the zoom particles"));

        pGrid[0]->gatherParticleList(zoomParticleArray);

        // find boundaries
        int x0, x1, y0, y1, z0, z1;
        int x,y,z;

        x0=y0=z0=n[0];
        x1=y1=z1=0;
        Grid<MyFloat> g(n[0]);

        // TO DO: wrap the box sensibly

        for(size_t i=0; i<zoomParticleArray.size(); i++) {
            g.getCoordinates(zoomParticleArray[i],x,y,z);
            if(x<x0) x0=x;
            if(y<y0) y0=y;
            if(z<z0) z0=z;
            if(x>x1) x1=x;
            if(y>y1) y1=y;
            if(z>z1) z1=z;
        }

        // Now see if the zoom the user chose is OK
        int n_user = n[0]/zoomfac;
        if((x1-x0)>n_user || (y1-y0)>n_user || (z1-z0)>n_user) {
            throw(std::runtime_error("Zoom particles do not fit in specified sub-box. Decrease zoom, or choose different particles. (NB wrapping not yet implemented)"));
        }
        // At this point we know things fit. All we need to do is choose
        // the correct offset to get the particles near the centre of the
        // zoom box.

        // Here is the bottom left of the box:
        x=(x0+x1)/2-n[0]/(2*zoomfac);
        y=(y0+y1)/2-n[0]/(2*zoomfac);
        z=(z0+z1)/2-n[0]/(2*zoomfac);

        if(x<0) x=0;
        if(y<0) y=0;
        if(z<0) z=0;
        if(x>n[0]-n[0]/zoomfac) x = n[0]-n[0]/zoomfac;
        if(y>n[0]-n[0]/zoomfac) y = n[0]-n[0]/zoomfac;
        if(z>n[0]-n[0]/zoomfac) z = n[0]-n[0]/zoomfac;


        x_off[1] = x_off[0] + x*dx[0];
        y_off[1] = y_off[0] + y*dx[0];
        z_off[1] = z_off[0] + z*dx[0];



        initZoom();

    }

    void initZoom() {
        zoomfac = (boxlen[0]/n[0])/(boxlen[1]/n[1]);
        dx[1] = boxlen[1]/n[1];
        cout << "Initialized a zoom region:" << endl;
        cout << "  Subbox length   = " << boxlen[1] << " Mpc/h" << endl;
        cout << "  n[1]            = " << n[1] << endl;
        cout << "  dx[1]           = " << dx[1] << endl;
        cout << "  Zoom factor     = " << zoomfac << endl;
        cout << "  Low-left corner = " << x_off[1] << ", " << y_off[1] << ", " << z_off[1] << endl;
        cout << "  Num particles   = " << zoomParticleArray.size() << endl;
        initGrid(1);

        updateParticleMapper();

        cout << "  Total particles = " << pMapper->size() << endl;
    }


    void setOutputMode(int in) {
        out = in; // the joy of writing this line is substantial

        if(out!=0 && out!=1 && out!=2)
            throw runtime_error("Wrong output format, choose 0 (HDF5), 1 (Gadget) or 2 (both)");

#ifndef HAVE_HDF5
        if(out!=1)
            throw runtime_error("Not compiled with HDF5. Only output=1 is allowed in this case!");
#endif

        updateParticleMapper();

    }

    void setSeed(int in) {
        randomFieldGenerator.seed(in);
    }

    void setSeedFourier(int in) {
        randomFieldGenerator.seed(in);
        randomFieldGenerator.setDrawInFourierSpace(true);
        randomFieldGenerator.setReverseRandomDrawOrder(false);
    }

    void setSeedFourierReverseOrder(int in) {
        randomFieldGenerator.seed(in);
        randomFieldGenerator.setDrawInFourierSpace(true);
        randomFieldGenerator.setReverseRandomDrawOrder(true);
    }

    void setExactPowerSpectrumEnforcement() {
        fieldManager.setExactPowerSpectrumEnforcement(true);
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
        updateParticleMapper();
        if(!prepared)
            prepare(); //compatibility with old paramfiles
    }

    string make_base( string basename, int level=0){
        ostringstream nult;
        if(inname.size()==0) {
            nult << basename<<"IC_iter_" << floatinfo<MyFloat>::name << "_z"<< cosmology.redshift <<"_"<<n[level]<<"_L" << boxlen[level];

        } else {
            if (level==0)
                nult << basename << "/" << inname;
            else
                nult << basename << "/" << inname << "_" << level;
        }
        return nult.str();
    }

    void readCamb() {
        spectrum.read(incamb);
        updateFieldManager();
    }

    void updateFieldManager() {
        if(spectrum.isUsable()) {
            for_each_level(level) {
                if(fieldManager.getNumLevels()<=level)
                    fieldManager.addLevel(spectrum.getPowerSpectrumForGrid(this->cosmology, *(this->pGrid[level])),
                                          this->pGrid[level],
                                          this->nPartLevel[level]);
            }
        }
    }


    virtual void drawRandom() {
        randomFieldGenerator.draw();
    }

    virtual void zeroLevel(int level) {
        cerr << "*** Warning: your script calls zeroLevel("<<level <<"). This is intended for testing purposes only!" << endl;
        fieldManager.zeroLevel(level);
    }

    template<typename T>
    void interpolateIntoLevel(const GridPtrType &pGrid_dest, const GridPtrType &pGrid_src,
                              T &pField_x_dest, T &pField_x_src) {
        pGrid_dest->addFieldFromDifferentGrid(*pGrid_src, pField_x_src, pField_x_dest);
    }

    virtual void interpolateIntoLevel(int level) {
        if(level<=0)
            throw std::runtime_error("Trying to interpolate onto the top-level grid");

        auto pGrid_l = pGrid[level];
        auto pGrid_p = pGrid[level-1];

        pGrid_l->addFieldFromDifferentGrid(*pGrid_p);
    }





    virtual void applyPowerSpec() {
        fieldManager.applyCovarianceToWhiteNoiseField();
    }

    template<typename T>
    void dumpGridData(int level, const T & data) {
        assert(data.size()==pGrid[level]->size3);
        ostringstream filename;
        filename << indir << "/grid-" << level << ".npy";

        const int dim[3] = { n[level],n[level],n[level] };
        numpy::SaveArrayAsNumpy(filename.str(), false, 3, dim, data.data());

        filename.str("");
        filename << indir << "/grid-info-" << level << ".txt";

        ofstream ifile;
        ifile.open(filename.str());
        cerr << "Writing to " << filename.str() << endl;

        ifile << x_off[level] << " " << y_off[level] << " " << z_off[level] << " " << boxlen[level] << endl;
        ifile << "The line above contains information about grid level " << level << endl;
        ifile << "It gives the x-offset, y-offset and z-offset of the low-left corner and also the box length" << endl;
        ifile.close();
    }

    virtual void saveTipsyArray(string fname) {
        saveFieldTipsyArray(fname, pMapper);
    }

    virtual void dumpGrid(int level=0) {
        dumpGridData(level, pGrid[level]->getFieldReal());
    }


    virtual void dumpPS(int level=0) {
        powsp_noJing(n[level], pGrid[level]->getFieldFourier(),
                     fieldManager.getCovariance(level),
                     (base+"_"+((char)(level+'0'))+".ps").c_str(), boxlen[level]);
    }

    virtual void zeldovichForLevel(int level) {
        //Zeldovich approx.

        MyFloat hfac=1.*100.*sqrt(cosmology.OmegaM0 / cosmology.scalefactor / cosmology.scalefactor / cosmology.scalefactor + cosmology.OmegaLambda0)*sqrt(
                cosmology.scalefactor);
        //this should be f*H(t)*a, but gadget wants vel/sqrt(a), so we use H(t)*sqrt(a)
        //TODO: hardcoded value of f=1 is inaccurate, but fomega currently gives wrong nults

        MyFloat pmass=27.78* cosmology.OmegaM0 *powf(boxlen[level]/(MyFloat)(n[level]),3.0);

        pGrid[level]->zeldovich(hfac,pmass);

    }



    virtual void zeldovich() {
        if(pGrid.size()==0) {
          throw std::runtime_error("Trying to apply zeldovich approximation, but no grids have been created");
        } else if(pGrid.size()==1) {
            zeldovichForLevel(0);
        } else {
            cerr << "Zeldovich approximation on successive levels...";
            zeldovichForLevel(0);
            zeldovichForLevel(1);
            cerr << "done." << endl;

            cerr << "Interpolating low-frequency information into zoom region...";
            interpolateIntoLevel(1);
            cerr << "done." << endl;

            cerr << "Re-introducing high-k modes into low-res region...";
            fieldManager.recombineLevel0();
            zeldovichForLevel(0);
            cerr << "done." << endl;

        }
    }

    void setInputMapper(std::string fname) {
        DummyIC<MyFloat> pseudoICs(this);
        auto dispatch = interpreter.specify_instance(pseudoICs);
        ifstream inf;
        inf.open(fname);


        if(!inf.is_open())
        throw std::runtime_error("Cannot open IC paramfile for relative_to command");
        cerr << "******** Running commands in" << fname << " to work out relationship ***********" << endl;

        ChangeCwdWhileInScope temporary(getDirectoryName(fname));

        dispatch.run_loop(inf);
        cerr << *(pseudoICs.pMapper) << endl;
        cerr << "******** Finished with" << fname << " ***********" << endl;
        pInputMapper = pseudoICs.pMapper;

    }

    std::shared_ptr<Grid<MyFloat>> getGridWithOutputOffset(int level=0) {
        auto gridForOutput = pGrid[level];

        if(xOffOutput!=0 || yOffOutput!=0 || zOffOutput!=0) {
            gridForOutput = std::make_shared<OffsetGrid<MyFloat>>(pGrid[0],xOffOutput,yOffOutput,zOffOutput);
        }
        return gridForOutput;
    }

    void updateParticleMapper() {

      if(pGrid.size()==0)
        return;

      // make a basic mapper for the base level grid
      pMapper = std::shared_ptr<ParticleMapper<MyFloat>>(new OneLevelParticleMapper<MyFloat>(
            getGridWithOutputOffset(0)
      ));


      if(pGrid.size()>=3) {
        // possible future enhancement, but for now...
        throw runtime_error("Don't know how to set up a mapper for more than one level of refinement");
      }

      if(pGrid.size()==2) {
        // it's a zoom!
        auto pMapperLevel1 = std::shared_ptr<ParticleMapper<MyFloat>>(new OneLevelParticleMapper<MyFloat>(getGridWithOutputOffset(1)));

        pMapper = std::shared_ptr<ParticleMapper<MyFloat>>(
            new TwoLevelParticleMapper<MyFloat>(pMapper, pMapperLevel1, zoomParticleArray, zoomfac*zoomfac*zoomfac));
      }

      if(cosmology.OmegaBaryons0 >0) {

          // Add gas only to the deepest level. Pass the whole pGrid
          // vector if you want to add gas to every level.
          auto gasMapper = pMapper->addGas(cosmology.OmegaBaryons0 / cosmology.OmegaM0,
                                          {pGrid.back()});

          bool gasFirst = gadgetformat==4;

          // graft the gas particles onto the start of the map
          if(gasFirst)
              pMapper = std::make_shared<AddGasMapper<MyFloat>>(
                      gasMapper.first, gasMapper.second, true);
          else
              pMapper = std::make_shared<AddGasMapper<MyFloat>>(
                      gasMapper.second, gasMapper.first, false);

      }

      // potentially resample the lowest-level DM grid. Again, this is theoretically
      // more flexible if you pass in other grid pointers.
      if(supersample>1)
          pMapper = pMapper->superOrSubSampleDM(supersample, {pGrid.back()},true);

      if(subsample>1)
        pMapper = pMapper->superOrSubSampleDM(subsample, {pGrid[0]}, false);

    }

    void clearAndDistributeParticleList() {

        if(pInputMapper!=nullptr) {
          pInputMapper->clearParticleList();
          pInputMapper->distributeParticleList(genericParticleArray);
        }
        else
        {
          pMapper->clearParticleList();
          pMapper->distributeParticleList(genericParticleArray);
        }
    }



    virtual void write() {

        zeldovich();

        cerr << "Write, ndm=" << pMapper->size_dm() << ", ngas=" << pMapper->size_gas() << endl;
        cerr << (*pMapper);

        if (gadgetformat==3 || gadgetformat==2)
            SaveGadget(base + ".gadget", boxlen[0], pMapper, cosmology, gadgetformat);

        else if (gadgetformat==4)
            SaveTipsy(base + ".tipsy", boxlen[0], pMapper, cosmology);

        else{ throw std::runtime_error("Invalid value for gadgetformat!");}
    }

    void makeInitialRealizationWithoutConstraints() {
        for_each_level(level) nPartLevel[level] = ((long)n[level]*n[level])*n[level];

        base = make_base(indir);

        readCamb();
        drawRandom();

        applyPowerSpec();
    }


    virtual void prepare() {
        if(prepared)
            throw(std::runtime_error("Called prepare, but grid is already prepared for constraints"));

        makeInitialRealizationWithoutConstraints();

        prepared=true;

    }



protected:

    int deepestLevelWithParticlesSelected() {
        if(pGrid.size()>1 && pGrid[1]->estimateParticleListSize()>0) return 1; else return 0;
    }

    int deepestLevel() {
        if(pGrid.size()>1) return 1; else return 0;
    }

    MyFloat get_wrapped_delta(MyFloat x0, MyFloat x1) {
        return pGrid[0]->getWrappedDelta(x0,x1);
    }


    void getCentre() {
        x0 = 0; y0 = 0; z0 =0;

        int level = deepestLevelWithParticlesSelected();

        MyFloat xa,ya,za, xb, yb, zb;

        std::vector<size_t> particleArray;
        pGrid[level]->gatherParticleList(particleArray);

        pGrid[level]->getCentroidLocation(particleArray[0],xa,ya,za);

        for(size_t i=0;i<particleArray.size();i++) {
            pGrid[level]->getCentroidLocation(particleArray[i],xb,yb,zb);
            x0+=get_wrapped_delta(xa,xb);
            y0+=get_wrapped_delta(ya,yb);
            z0+=get_wrapped_delta(za,zb);
        }
        x0/=particleArray.size();
        y0/=particleArray.size();
        z0/=particleArray.size();
        x0+=xa;
        y0+=ya;
        z0+=za;
        cerr << "Centre of region is " << x0 << " " << y0 << " " << z0 << endl;
    }




    void appendParticleIdFile(std::string filename) {

        cerr << "Loading " << filename << endl;

        getBuffer(genericParticleArray, filename);

        cerr << "  -> total number of particles is " << genericParticleArray.size() << endl;

        clearAndDistributeParticleList();
    }

    void loadParticleIdFile(std::string filename) {
        genericParticleArray.clear();
        appendParticleIdFile(filename);
    }


    vector<complex<MyFloat>> calcConstraint(string name_in, bool kspace = true) {
        return constraintGenerator.calcConstraintForAllLevels(name_in, kspace);
    }



public:


    void loadID(string fname) {
        appendParticleIdFile(fname);
        getCentre();
    }

    void appendID(string fname) {
        appendParticleIdFile(fname);
        getCentre();
    }

    void dumpID(string fname) {
        std::vector<size_t> results;
        pMapper->gatherParticleList(results);
        dumpBuffer(results, fname);
    }

    void centreParticle(long id) {
        pGrid[0]->getCentroidLocation(id,x0,y0,z0);
    }

    void selectNearest() {
        auto grid = pGrid[deepestLevel()];
        pMapper->clearParticleList();
        size_t id = grid->getClosestIdNoWrap(x0,y0,z0);
        grid->distributeParticleList({id});

    }

    void selectSphere(float radius) {
        MyFloat r2 = radius*radius;
        MyFloat delta_x, delta_y, delta_z, r2_i;
        MyFloat xp,yp,zp;

        genericParticleArray.clear();

        for_each_level(level) {
            std::vector<size_t> particleArray;

            for(size_t i=0;i<this->nPartLevel[level];i++) {
                pGrid[level]->getCentroidLocation(i,xp,yp,zp);
                delta_x = get_wrapped_delta(xp,x0);
                delta_y = get_wrapped_delta(yp,y0);
                delta_z = get_wrapped_delta(zp,z0);
                r2_i = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
                if(r2_i<r2)
                  particleArray.push_back(i);
            }

            pGrid[level]->clearParticleList();
            pGrid[level]->distributeParticleList(particleArray);

        }

    }


    void setCentre(MyFloat xin, MyFloat yin, MyFloat zin) {
            x0=xin;
            y0=yin;
            z0=zin;
    }

    void reorderBuffer() {

        cout << "Reordering buffer radially..." << endl;
        cout << " [taking centre = " << x0 << " " << y0 << " " <<z0 << "]" << endl;

        std::vector<MyFloat> r2(genericParticleArray.size());
        MyFloat delta_x, delta_y, delta_z;

        std::vector<size_t> index(genericParticleArray.size());

        for(size_t i=0;i<genericParticleArray.size();i++) {
            pGrid[0]->getCentroidLocation(i,delta_x, delta_y, delta_z);
            delta_x = get_wrapped_delta(delta_x,x0);
            delta_y = get_wrapped_delta(delta_y,y0);
            delta_z = get_wrapped_delta(delta_z,z0);
            r2[i] = delta_x*delta_x+delta_y*delta_y+delta_z*delta_z;
            index[i]=i;
        }

        // Now sort the index array
        std::sort(index.begin(),index.end(),
             [&r2](size_t i1, size_t i2) { return r2[i1]<r2[i2]; } );

        // Turn the index array into something pointing to the particles
        for(size_t i=0; i<genericParticleArray.size(); i++) {
            index[i] = genericParticleArray[index[i]];
        }

        // Copy back into the particle array
        for(size_t i=0; i<genericParticleArray.size(); i++) {
            genericParticleArray[i] = index[i];
        }

    }

    void truncateBuffer(float x) {
        if(x<0 ) {
            cerr << "Truncate command takes a fraction between 0 and 1 or a number>=1" << endl;
            exit(0);
        }
        if(x<1)
            genericParticleArray.resize(((int)(genericParticleArray.size()*x)));
        else
            genericParticleArray.resize(((int)x));
    }

    void calculate(string name) {
        auto vecs = calcConstraint(name);
        float val = float((fieldManager.v1_dot_y(vecs)).real());
        cout << name << ": calculated value = " <<  val << endl;
    }

    virtual void constrain(string name, string type, float value) {
        bool relative=false;
        if (strcasecmp(type.c_str(),"relative")==0) {
            relative=true;
        } else if (strcasecmp(type.c_str(),"absolute")!=0) {
            throw runtime_error("Constraint type must be either 'relative' or 'absolute'");
        }

        std::complex<MyFloat> constraint = value;
        auto vec = calcConstraint(name);
        std::complex<MyFloat> initv = fieldManager.v1_dot_y(vec);

        if(relative) constraint*=initv;

        cout << name << ": initial value = " << initv << ", constraining to " << constraint << endl;
        constraintApplicator.add_constraint(std::move(vec), constraint, initv);

    }

    void cov() {
        constraintApplicator.print_covariance();
    }


    virtual void fixConstraints() {
        constraintApplicator.prepare();
        constraintApplicator.applyConstraints();
    }

    virtual void done() {
        cerr << "BEFORE constraints: chi^2=" << fieldManager.get_field_chi2() << endl;
        fixConstraints();
        cerr << "AFTER  constraints: chi^2=" << fieldManager.get_field_chi2() << endl;
        write();
    }

    void reverse() {
        for_each_level(level) {
            auto & field = pGrid[level]->getField();
            for(size_t i=0; i<this->nPartLevel[level]; i++)
                field[i]=-field[i];
        }
    }

    void reseedSmallK(MyFloat kmax, int seed) {

        MyFloat k2max = kmax*kmax;

        // take a copy of all the fields
        std::vector<FieldType> fieldCopies;
        for_each_level(level) fieldCopies.emplace_back(pGrid[level]->getFieldFourier());

        // remake the fields with the new seed
        randomFieldGenerator.seed(seed);
        makeInitialRealizationWithoutConstraints();

        // copy back the old field
        for_each_level(level) {
            auto & fieldOriginal = fieldCopies[level];
            auto & field = pGrid[level]->getFieldFourier();
            const auto & grid = *(this->pGrid[level]);
            MyFloat k2;
            for(size_t i=0; i<this->nPartLevel[level]; i++) {
                k2 = grid.getKSquared(i);
                if(k2>k2max && k2!=0) {
                    field[i]=fieldOriginal[i];
                }
            }
        }

    }

    void reverseSmallK(MyFloat kmax) {

        MyFloat k2max = kmax*kmax;


        for_each_level(level) {
            MyFloat k2_g_min=std::numeric_limits<MyFloat>::max();
            MyFloat k2_g_max=0.0;
            size_t modes_reversed=0;
            size_t tot_modes =pGrid[level]->size3;
            auto & field = pGrid[level]->getFieldFourier();
            const auto & grid = *(this->pGrid[level]);
            MyFloat k2;
            for(size_t i=0; i<this->nPartLevel[level]; i++) {
                k2 = grid.getKSquared(i);
                if(k2<k2max && k2!=0) {
                    field[i]=-field[i];
                    modes_reversed++;
                }
                if(k2<k2_g_min && k2!=0)
                    k2_g_min = k2;
                if(k2>k2_g_max)
                    k2_g_max = k2;
            }
            cerr << "reverseSmallK: k reversal at " << sqrt(k2max) << "; grid was in range " << sqrt(k2_g_min) << " to " << sqrt (k2_g_max) << endl;
            cerr << "               modes reversed = " << modes_reversed << " of " << tot_modes << endl;
        }

    }



};
#endif
