#ifndef _IO_INCLUDED
#define _IO_INCLUDED

#include <cassert>
#include <vector>

#include "mapper.hpp"

#ifdef HAVE_HDF5
hid_t hdf_float = H5Tcopy (H5T_NATIVE_DOUBLE);
hid_t hdf_double = H5Tcopy (H5T_NATIVE_DOUBLE);
#endif

// #endif



using namespace std;



size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream) //stolen from Gadget
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) has occured.\n");
      fflush(stdout);
    }
  return nwritten;
}

struct io_header_tipsy
{
    double scalefactor;
    int n;
    int ndim;
    int ngas;
    int ndark;
    int nstar;
} header_tipsy;

struct io_tipsy_dark
{
    float mass,x,y,z,vx,vy,vz,eps,phi;
};

struct io_tipsy_gas
{
    float mass,x,y,z,vx,vy,vz,rho,temp,eps,metals,phi;
};



struct io_header_2 //header for gadget2
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      nPartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the nPartTotal[2] stores the result of a division of the particle number by 2^32, while nPartTotal[1] holds the remainder. */
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header2;


struct io_header_3 //header for gadget3
{
  int npart[6];			/*!< number of particles of each type in this file */
  double mass[6];		/*!< mass of particles of each type. If 0, then the masses are explicitly
                   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;			/*!< time of snapshot file */
  double redshift;		/*!< redshift of snapshot file */
  int flag_sfr;			/*!< flags whether the simulation was including star formation */
  int flag_feedback;		/*!< flags whether feedback was included (obsolete) */
  unsigned int nPartTotal[6];	/*!< total number of particles of each type in this snapshot. This can be
                   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;		/*!< flags whether cooling was included  */
  int num_files;		/*!< number of files in multi-file snapshot */
  double BoxSize;		/*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;		/*!< matter density in units of critical density */
  double OmegaLambda;		/*!< cosmological constant parameter */
  double HubbleParam;		/*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;		/*!< flags whether the file contains formation times of star particles */
  int flag_metals;		/*!< flags whether the file contains metallicity values for gas and star
                   particles */
  unsigned int nPartTotalHighWord[6];	/*!< High word of the total number of particles of each type (see header2)*/
  int flag_entropy_instead_u;	/*!< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;	/*!< flags that snapshot contains double-precision instead of single precision */

  int flag_ic_info;             /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
                                     or whether they ocontains 2nd order lagrangian perturbation theory initial conditions.
                                     For snapshots files, the value informs whether the simulation was evolved from
                                     Zeldoch or 2lpt ICs. Encoding is as follows:
                                        FLAG_ZELDOVICH_ICS     (1)   - IC file based on Zeldovich
                                        FLAG_SECOND_ORDER_ICS  (2)   - Special IC-file containing 2lpt masses
                                        FLAG_EVOLVED_ZELDOVICH (3)   - snapshot evolved from Zeldovich ICs
                                        FLAG_EVOLVED_2LPT      (4)   - snapshot evolved from 2lpt ICs
                                        FLAG_NORMALICS_2LPT    (5)   - standard gadget file format with 2lpt ICs
                                     All other values, including 0 are interpreted as "don't know" for backwards compatability.
                                 */
  float lpt_scalingfactor;      /*!< scaling factor for 2lpt initial conditions */

  char fill[48];		/*!< fills to 256 Bytes */

} header3;				/*!< holds header for snapshot files */

#ifdef HAVE_HDF5
template<typename MyFloat>
void SaveHDF(const char* Filename, int n, io_header_3 header, MyFloat *p_Data1,  MyFloat *p_Data2, MyFloat *p_Data3, MyFloat *v_Data1,  MyFloat *v_Data2, MyFloat *v_Data3, const char *name1, const char *name2, const char *name3) //saves 3 arrays as datasets "name1", "name2" and "name3" in one HDF5 file
{
  hid_t       file_id, gidh, gid, dset_id,ndset_id,iddset_id, hdset_id,mt_id;   /* file and dataset identifiers */
  hid_t       nfilespace,dataspace_id, memspace, hboxfs,mt,idfilespace,IDmemspace,iddataspace_id;/* file and memory dataspace identifiers */
  hsize_t ncount[2],ncount2[2], noffset[2], h[1], stride[2],  idcount[2],idoffset[2],idstride[2],mta[1],idblock[2];
  h[0]=1;
  int m_nGrid[3]={n,n,n};
  mta[0]=6;

  file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  gidh=H5Gcreate(file_id, "Header", 100 ); //100 gives max number of arguments (?)

   MyFloat *box=new MyFloat[1];
   MyFloat *mtt=new MyFloat[6];
   int *boxi=new int[1];
   mtt[0]=header.mass[0];
   mtt[1]=header.mass[1];
   mtt[2]=header.mass[2];
   mtt[3]=header.mass[3];
   mtt[4]=header.mass[4];
   mtt[5]=header.mass[5];
   mt = H5Screate_simple( 1, mta, NULL );

   mt_id = H5Acreate( gidh, "MassTable", hdf_float, mt, H5P_DEFAULT );
   H5Awrite( mt_id, hdf_float, mtt);
   H5Aclose(mt_id);

   hboxfs = H5Screate_simple( 1, h, NULL );
   box[0]=header.Omega0;
   hdset_id = H5Acreate( gidh, "OM", hdf_float, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, hdf_float, box);
   H5Aclose(hdset_id);

   box[0]=header.OmegaLambda;
   hdset_id = H5Acreate( gidh, "OLambda", hdf_float, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, hdf_float, box);
   H5Aclose(hdset_id);

   box[0]=header.BoxSize;
   hdset_id = H5Acreate( gidh, "BoxSize", hdf_float, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, hdf_float, box);
   H5Aclose(hdset_id);

   box[0]=header.redshift;
   hdset_id = H5Acreate( gidh, "Redshift", hdf_float, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, hdf_float, box);
   H5Aclose(hdset_id);

   box[0]=header.time;
   hdset_id = H5Acreate( gidh, "Time", hdf_float, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, hdf_float, box);
   H5Aclose(hdset_id);

   //box[0]=0.; //fix: read as input
   //hdset_id = H5Acreate( gidh, "fnl", hdf_float, hboxfs, H5P_DEFAULT );
   //H5Awrite( hdset_id, hdf_float, box);
   //H5Aclose(hdset_id);

   box[0]=0.817;//fix: read as input
   hdset_id = H5Acreate( gidh, "sigma8", hdf_float, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, hdf_float, box);
   H5Aclose(hdset_id);

   long *boxu=(long*)calloc(6,sizeof(long));
   boxu[1]=(long)(m_nGrid[0]*m_nGrid[1]*m_nGrid[2]);
   hdset_id = H5Acreate( gidh, "NumPart_Total", H5T_NATIVE_UINT, mt, H5P_DEFAULT );
   H5Awrite( hdset_id, H5T_NATIVE_UINT, boxu);
   H5Aclose(hdset_id);

//    hdset_id = H5Acreate( gidh, "NumPart_ThisFile", H5T_NATIVE_UINT, mt, H5P_DEFAULT );
//    H5Awrite( hdset_id, H5T_NATIVE_UINT, boxu);
//    H5Aclose(hdset_id);

   boxu[1]=header.nPartTotalHighWord[1];
   hdset_id = H5Acreate( gidh, "NumPart_Total_HighWord", H5T_NATIVE_UINT, mt, H5P_DEFAULT );
   H5Awrite( hdset_id, H5T_NATIVE_UINT, boxu);
   H5Aclose(hdset_id);

   boxi[0]=header.flag_ic_info;
   hdset_id = H5Acreate( gidh, "Flag_IC_Info", H5T_NATIVE_INT, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, H5T_NATIVE_INT, boxi);
   H5Aclose(hdset_id);

   boxi[0]=header.flag_doubleprecision;
   hdset_id = H5Acreate( gidh, "flag_doubleprecision", H5T_NATIVE_INT, hboxfs, H5P_DEFAULT );
   H5Awrite( hdset_id, H5T_NATIVE_INT, boxi);
   H5Aclose(hdset_id);

//    boxi[0]=header.num_files;
//    hdset_id = H5Acreate( gidh, "NumFilesPerSnapshot", H5T_NATIVE_INT, hboxfs, H5P_DEFAULT );
//    H5Awrite( hdset_id, H5T_NATIVE_INT, boxi);
//    H5Aclose(hdset_id);


   H5Sclose(mt);
   H5Sclose(hboxfs);

    //close "header" group
       H5Gclose(gidh);

  gid=H5Gcreate(file_id, "PartType1", 100 );
  idcount[0]=m_nGrid[0]*m_nGrid[1]*m_nGrid[2];
  idcount[1]=1;
  idfilespace = H5Screate_simple( 2, idcount, NULL );
  iddset_id = H5Dcreate( gid, name3, H5T_NATIVE_LONG, idfilespace, H5P_DEFAULT );
  H5Sclose(idfilespace);

  IDmemspace = H5Screate_simple( 2, idcount, NULL );

  idoffset[0]=0;
  idoffset[1]=0;
  idstride[0]=1;
  idstride[1]=1;
  idblock[0]=1;
  idblock[1]=1;

  long i;
  long* ID=(long*)calloc(n*n*n,sizeof(long));
  for (i=0; i< (long)(n*n*n); i++){ID[i]=i;}

  iddataspace_id= H5Dget_space(iddset_id);
  H5Sselect_hyperslab(iddataspace_id, H5S_SELECT_SET, idoffset, idstride, idcount, idblock);
  H5Dwrite( iddset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ID );

  free(ID);

  H5Sclose(iddataspace_id);
  H5Sclose(IDmemspace);
  H5Dclose(iddset_id);

  ncount[0]=m_nGrid[0]*m_nGrid[1]*m_nGrid[2]; //total number of particles
  ncount[1]=3; //3 components of position vector
  nfilespace = H5Screate_simple( 2, ncount, NULL ); //works! for (n**3, (x1,x2,x3) ) array (part pos./vel.)
  dset_id = H5Dcreate( gid, name1, hdf_float, nfilespace, H5P_DEFAULT );
  ndset_id = H5Dcreate( gid, name2, hdf_float, nfilespace, H5P_DEFAULT );
  H5Sclose(nfilespace);

  noffset[0] = 0;
  noffset[1] = 0;

  ncount2[0]=m_nGrid[0]*m_nGrid[1]*m_nGrid[2];
  ncount2[1]=1;

  memspace = H5Screate_simple( 2, ncount2, NULL );

  dataspace_id= H5Dget_space(dset_id);

  stride[0]=1;
  stride[1]=1;

  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, noffset, stride, ncount2, NULL);
    H5Dwrite( dset_id, hdf_float, memspace, dataspace_id, H5P_DEFAULT, p_Data1 );

    noffset[0] = 0;
    noffset[1] = 1;

    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, noffset, stride, ncount2, NULL);
    H5Dwrite( dset_id, hdf_float, memspace, dataspace_id, H5P_DEFAULT, p_Data2 );

    noffset[0] = 0;
    noffset[1] = 2;

    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, noffset, stride, ncount2, NULL);
    H5Dwrite( dset_id, hdf_float, memspace, dataspace_id, H5P_DEFAULT, p_Data3 );

    noffset[0] = 0;
    noffset[1] = 0;
   H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, noffset, NULL, ncount2, NULL);
   H5Dwrite( ndset_id, hdf_float, memspace, dataspace_id, H5P_DEFAULT, v_Data1 );

   noffset[0] = 0;
    noffset[1] = 1;
   H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, noffset, NULL, ncount2, NULL);
   H5Dwrite( ndset_id, hdf_float, memspace, dataspace_id, H5P_DEFAULT, v_Data2 );

   noffset[0] = 0;
    noffset[1] = 2;
   H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, noffset, NULL, ncount2, NULL);
   H5Dwrite( ndset_id, hdf_float, memspace, dataspace_id, H5P_DEFAULT, v_Data3 );

  H5Sclose(memspace);
  H5Dclose(dset_id);
  H5Dclose(ndset_id);
  H5Sclose(dataspace_id);

  //close "particle" group
    H5Gclose(gid);

      //close overall file
         H5Fclose(file_id);

    free(box);
    free(boxi);
    free(boxu);
    free(mtt);

}

template<typename MyFloat>
int save_phases(complex<MyFloat> *phk, MyFloat* ph, complex<MyFloat> *delta, int n, const char *name){

      MyFloat *helper=(MyFloat*)calloc(n*n*n, sizeof(MyFloat));
      int i;
      for(i=0;i<n*n*n;i++){helper[i]=delta[i].real();}

    hid_t     idfilespace,  file_id, gid, iddset_id,IDmemspace,iddataspace_id ;
    hsize_t idcount[2];
    file_id = H5Fcreate( name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    gid=H5Gcreate(file_id, "Potential", 100 );
    idcount[0]=n*n*n;
    idcount[1]=2;
    idfilespace = H5Screate_simple( 2, idcount, NULL );
    iddset_id = H5Dcreate( gid, "Psi_k", hdf_float, idfilespace, H5P_DEFAULT );
    H5Sclose(idfilespace);

    IDmemspace = H5Screate_simple( 2, idcount, NULL );

    iddataspace_id= H5Dget_space(iddset_id);
    H5Dwrite( iddset_id, hdf_float, H5S_ALL, H5S_ALL, H5P_DEFAULT, phk );

    H5Sclose(iddataspace_id);
    H5Sclose(IDmemspace);
    H5Dclose(iddset_id);

    idcount[1]=1;
    idfilespace = H5Screate_simple( 2, idcount, NULL );
    iddset_id = H5Dcreate( gid, "Psi_r", hdf_float, idfilespace, H5P_DEFAULT );
    H5Sclose(idfilespace);

    IDmemspace = H5Screate_simple( 2, idcount, NULL );

    iddataspace_id= H5Dget_space(iddset_id);
    H5Dwrite( iddset_id, hdf_float, H5S_ALL, H5S_ALL, H5P_DEFAULT, ph );

    H5Sclose(iddataspace_id);
    H5Sclose(IDmemspace);
    H5Dclose(iddset_id);

    idfilespace = H5Screate_simple( 2, idcount, NULL );
    iddset_id = H5Dcreate( gid, "delta_r", hdf_float, idfilespace, H5P_DEFAULT );
    H5Sclose(idfilespace);

    IDmemspace = H5Screate_simple( 2, idcount, NULL );

    iddataspace_id= H5Dget_space(iddset_id);
    H5Dwrite( iddset_id, hdf_float, H5S_ALL, H5S_ALL, H5P_DEFAULT, helper );

    H5Sclose(iddataspace_id);
    H5Sclose(IDmemspace);
    H5Dclose(iddset_id);

    //close "particle" group
    H5Gclose(gid);

        //close overall file
        H5Fclose(file_id);

    return 0;
}
#endif

template<typename MyFloat>
void SaveGadget2(const char *filename, long nPart, io_header_2 header1, MyFloat* Pos1, MyFloat* Vel1, MyFloat* Pos2, MyFloat* Vel2, MyFloat* Pos3, MyFloat* Vel3, MyFloat *Mass=NULL) {
  FILE* fd = fopen(filename, "w");
  if(!fd) throw std::runtime_error("Unable to open file for writing");
  MyFloat* Pos=(MyFloat*)calloc(3,sizeof(MyFloat));
  int dummy;
  long i;
  //header block
  dummy= sizeof(header1);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  my_fwrite(&header1, sizeof(header1), 1, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  //position block
  dummy=sizeof(MyFloat)*(long)(nPart)*3; //this will be 0 or some strange number for n>563; BUT: gagdget does not actually use this value; it gets the number of particles from the header
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i=0;i<nPart;i++){
    Pos[0]=Pos1[i];
    Pos[1]=Pos2[i];
    Pos[2]=Pos3[i];
    my_fwrite(Pos,sizeof(MyFloat),3,fd);
  }
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  //velocity block
  //dummy=sizeof(MyFloat)*3*nPart; //this will be 0 or some strange number for n>563; BUT: gagdget does not actually use this value; it gets the number of particles from the header
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i=0;i<nPart;i++){
    Pos[0]=Vel1[i];
    Pos[1]=Vel2[i];
    Pos[2]=Vel3[i];
   my_fwrite(Pos,sizeof(MyFloat),3,fd);
  }

  free(Pos);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  dummy = sizeof(long) * nPart; //here: gadget just checks if the IDs are ints or long longs; still the number of particles is read from the header file
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i=0;i<nPart;i++){
  my_fwrite(&i, sizeof(long), 1, fd);

  }
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  if(Mass!=NULL) {
      dummy = sizeof(MyFloat)*nPart;
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      my_fwrite(Mass,sizeof(MyFloat),  nPart, fd);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
  }


  fclose(fd);

}

template<typename MyFloat>
void SaveGadget3(const char *filename, long n, io_header_3 header1, MyFloat* Pos1, MyFloat* Vel1, MyFloat* Pos2, MyFloat* Vel2, MyFloat* Pos3, MyFloat* Vel3, MyFloat *Mass=NULL) {

    FILE* fd = fopen(filename, "w");
    if(!fd) throw std::runtime_error("Unable to open file for writing");
    MyFloat* Pos=(MyFloat*)calloc(3,sizeof(MyFloat));
    int dummy;
    long i;

    //header block
    dummy= sizeof(header1);
    my_fwrite(&dummy, sizeof(dummy), 1, fd);
    my_fwrite(&header1, sizeof(header1), 1, fd);
    my_fwrite(&dummy, sizeof(dummy), 1, fd);

    //position block
    dummy=sizeof(MyFloat)*(long)(n)*3; //this will be 0 or some strange number for n>563; BUT: gagdget does not actually use this value; it gets the number of particles from the header
    my_fwrite(&dummy, sizeof(dummy), 1, fd);
    for(i=0;i<n;i++){
    Pos[0]=Pos1[i];
        Pos[1]=Pos2[i];
        Pos[2]=Pos3[i];
        my_fwrite(Pos,sizeof(MyFloat),3,fd);
     }

     my_fwrite(&dummy, sizeof(dummy), 1, fd);

    //velocity block
    my_fwrite(&dummy, sizeof(dummy), 1, fd);
    for(i=0;i<n;i++){
        Pos[0]=Vel1[i];
        Pos[1]=Vel2[i];
        Pos[2]=Vel3[i];
        my_fwrite(Pos,sizeof(MyFloat),3,fd);
       }

    my_fwrite(&dummy, sizeof(dummy), 1, fd);

    //particle block
    //long long ido;
    dummy = sizeof(long) * n; //here: gadget just checks if the IDs are ints or long longs; still the number of particles is read from the header file
    my_fwrite(&dummy, sizeof(dummy), 1, fd);
    for(i=0;i<n;i++){
         my_fwrite(&i, sizeof(long), 1, fd);
    }
    my_fwrite(&dummy, sizeof(dummy), 1, fd);


    if(Mass!=NULL) {
        dummy = sizeof(MyFloat)*n;
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
        my_fwrite(Mass,sizeof(MyFloat),  n, fd);
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
    }

    fclose(fd);
    free(Pos);

}

template<typename MyFloat>
void SaveTipsy(const std::string & filename,
               double Boxlength, double Om0, double Ol0, double hubble, double ain,
               shared_ptr<ParticleMapper<MyFloat>> pMapper) {

    // originally:
    // pmass in 1e10 h^-1 Msol
    // pos in Mpc h^-1
    // vel in km s^-1 a^1/2



    ofstream photogenic_file;

    double pos_factor  = 1./Boxlength;  // boxsize = 1
    double vel_factor  = ain/(sqrt(3./(8.*M_PI))*100*Boxlength);
    double mass_factor = 0.0; // calculated below shortly

    double min_mass=1.0/0;
    double max_mass=0.0;

    long n_gas = pMapper->size_gas();

    size_t iord=0;

    MyFloat x,y,z,vx,vy,vz,mass,tot_mass=0.0,eps;

    for(auto i=pMapper->begin(); i!=pMapper->end(); ++i) {
        i.getParticle(x,y,z,vx,vy,vz,mass,eps);
        if(min_mass>mass) min_mass=mass;
        if(max_mass<mass) max_mass=mass;
        tot_mass+=mass;
    }

    if(min_mass!=max_mass) {
        photogenic_file.open("photogenic.txt");
    }

    mass_factor = Om0/tot_mass; // tipsy convention: sum(mass)=Om0


    io_header_tipsy header;

    header.scalefactor = ain;
    header.n = pMapper->size();
    header.ndim = 3;
    header.ngas = pMapper->size_gas();
    header.ndark = pMapper->size_dm();
    header.nstar = 0;


    cout << "TIPSY parameters:" << endl;

    double dKpcUnit  = Boxlength*1000/hubble;
    double dMsolUnit = 1e10/hubble/mass_factor;
    double dKmsUnit  = sqrt(4.30211349e-6*dMsolUnit/(dKpcUnit));

    cout << "dKpcUnit: " <<  dKpcUnit << endl;
    cout << "dMsolUnit: " << dMsolUnit  << endl;
    cout << "hubble0: " << 0.1*hubble * dKpcUnit / dKmsUnit << endl;

    io_tipsy_dark dp;
    io_tipsy_gas gp;

    FILE* fd = fopen(filename.c_str(), "w");
    if(!fd) throw std::runtime_error("Unable to open file for writing");

    dp.phi = 0.0;
    gp.temp = 2.73/ain;
    gp.metals = 0.0;
    gp.rho = 0.0;


    fwrite(&header, sizeof(io_header_tipsy), 1, fd);


    for(auto i=pMapper->beginGas(); i!=pMapper->endGas(); ++i) {
        i.getParticle(x,y,z,vx,vy,vz,mass,eps);

        gp.x=x*pos_factor-0.5;
        gp.y=y*pos_factor-0.5;
        gp.z=z*pos_factor-0.5;

        gp.eps=eps*pos_factor;
        gp.vx=vx*vel_factor;
        gp.vy=vy*vel_factor;
        gp.vz=vz*vel_factor;
        gp.mass = mass*mass_factor;

        fwrite(&gp, sizeof(io_tipsy_gas), 1, fd);

    }

    for(auto i=pMapper->beginDm(); i!=pMapper->endDm(); ++i) {
        i.getParticle(x,y,z,vx,vy,vz,mass,eps);

        dp.x=x*pos_factor-0.5;
        dp.y=y*pos_factor-0.5;
        dp.z=z*pos_factor-0.5;

        dp.eps=eps*pos_factor;
        dp.vx=vx*vel_factor;
        dp.vy=vy*vel_factor;
        dp.vz=vz*vel_factor;
        dp.mass = mass*mass_factor;

        fwrite(&dp, sizeof(io_tipsy_dark), 1, fd);

    }
    fclose(fd);


}


template<typename MyFloat>
io_header_2 CreateGadget2Header(long nPartTotal, double pmass, double ain, double zin, double Boxlength, double Om0, double Ol0, double hubble)
{
    io_header_2 header2;
    header2.npart[0]=0;
    header2.npart[1]=nPartTotal;
    header2.npart[2]=0;
    header2.npart[3]=0;
    header2.npart[4]=0;
    header2.npart[5]=0;
    header2.mass[0]=0;
    header2.mass[1]=pmass;
    header2.mass[2]=0;
    header2.mass[3]=0;
    header2.mass[4]=0;
    header2.mass[5]=0;
    header2.time=ain;
    header2.redshift=zin;
    header2.flag_sfr=0;
    header2.flag_feedback=0;
    header2.nPartTotal[0]=0;
    header2.nPartTotal[1]=nPartTotal;
    header2.nPartTotal[2]=0;
    header2.nPartTotal[3]=0;
    header2.nPartTotal[4]=0;
    header2.nPartTotal[5]=0;
    header2.flag_cooling=0;
    header2.num_files=1;
    header2.BoxSize=Boxlength;
    header2.Omega0=Om0;
    header2.OmegaLambda=Ol0;
    header2.HubbleParam=hubble;
    return header2;
}

template<typename MyFloat>
io_header_3 CreateGadget3Header(long nPartTotal, double pmass, double ain, double zin, double Boxlength, double Om0, double Ol0, double hubble)
{
    io_header_3 header3;
    header3.npart[0]=0;
    header3.npart[1]=(unsigned int)(nPartTotal);
    header3.npart[2]=0;
    header3.npart[3]=0;
    header3.npart[4]=0;
    header3.npart[5]=0;
    header3.mass[0]=0;
    header3.mass[1]=pmass;
    header3.mass[2]=0;
    header3.mass[3]=0;
    header3.mass[4]=0;
    header3.mass[5]=0;
    header3.time=ain;
    header3.redshift=zin;
    header3.flag_sfr=0;
    header3.flag_feedback=0;
    header3.nPartTotal[0]=0;
    header3.nPartTotal[1]=(unsigned int)(nPartTotal);
    header3.nPartTotal[2]=0;
    header3.nPartTotal[3]=0;
    header3.nPartTotal[4]=0;
    header3.nPartTotal[5]=0;
    header3.flag_cooling=0;
    header3.num_files=1;
    header3.BoxSize=Boxlength;
    header3.Omega0=Om0;
    header3.OmegaLambda=Ol0;
    header3.HubbleParam=hubble;
    header3.flag_stellarage=0;	/*!< flags whether the file contains formation times of star particles */
    header3.flag_metals=0;		/*!< flags whether the file contains metallicity values for gas and star  particles */
    header3.nPartTotalHighWord[0]=0;
    header3.nPartTotalHighWord[1]=(unsigned int) (nPartTotal >> 32); //copied from Gadget3
    header3.nPartTotalHighWord[2]=0;
    header3.nPartTotalHighWord[3]=0;
    header3.nPartTotalHighWord[4]=0;
    header3.nPartTotalHighWord[5]=0;
    header3.flag_entropy_instead_u=0;	/*!< flags that IC-file contains entropy instead of u */
    header3.flag_doubleprecision=floatinfo<MyFloat>::doubleprecision;
    header3.flag_ic_info=1;
    return header3;
}


/*
INPUT ROUTINES
*/

template<typename T>
void getBuffer(std::vector<T> &store, std::string filename) {
    std::ifstream f(filename);
    if(!f.is_open())
        throw std::runtime_error("File "+filename+" not found");
    string line;
    while(!f.eof()) {
        T temp;
        if(f >> temp)
            store.push_back(temp);
    }
}

#endif // _IO_INCLUDED
