//
// Created by Andrew Pontzen on 18/11/2016.
//

#ifndef IC_GADGET_HPP
#define IC_GADGET_HPP

#include "src/io.hpp"

namespace io {
  namespace gadget {
    using std::cerr;
    using std::endl;

    template<typename GridDataType, typename FloatType=strip_complex<GridDataType>>
    void getParticleInfo(const particle::AbstractMultiLevelParticleGenerator<GridDataType> &generator,
                         particle::mapper::ParticleMapper<GridDataType> &mapper,
                         FloatType &min_mass, FloatType &max_mass,
                         FloatType &tot_mass, FloatType &gas_mass, size_t &ngas, size_t &nlow, size_t &nhigh) {

      min_mass = std::numeric_limits<FloatType>::max();
      max_mass = 0;
      gas_mass = 0;
      tot_mass = 0;
      ngas = 0;
      nlow = 0;
      nhigh = 0;

      FloatType mass;
      for (auto i = mapper.beginDm(generator); i != mapper.endDm(generator); ++i) {
        // progress("Pre-write scan file",iord, totlen);
        mass = i.getMass(); // sometimes can be MUCH faster than getParticle
        if (min_mass > mass) min_mass = mass;
        if (max_mass < mass) max_mass = mass;
        tot_mass += mass;
      }

      // end_progress();

      ngas = mapper.size_gas();
      if (ngas > 0) gas_mass = mapper.beginGas(generator).getMass();

      cerr << "gas mass and number of particles in info " << gas_mass << " " << ngas << endl;


      for (auto i = mapper.beginDm(generator); i != mapper.endDm(generator); ++i) {
        // progress("Pre-write scan file",iord, totlen);
        mass = i.getMass(); // sometimes can be MUCH faster than getParticle
        if (mass == min_mass) nhigh += 1;
        else if (mass == max_mass) nlow += 1;
        else {
          cerr << "else in mass " << min_mass << " " << max_mass << " " << mass << endl;
          continue;
        }

      }
      // end_progress();
    }

    size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream) //stolen from Gadget
    {
      size_t nwritten;

      if ((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb) {
        printf("I/O error (fwrite) has occured.\n");
        fflush(stdout);
      }
      return nwritten;
    }

    struct io_header_2 //header for gadget2
    {
      int npart[6];
      double mass[6];
      double time;
      double redshift;
      int flag_sfr;
      int flag_feedback;
      int nPartTotal[6];
      /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the nPartTotal[2] stores the result of a division of the particle number by 2^32, while nPartTotal[1] holds the remainder. */
      int flag_cooling;
      int num_files;
      double BoxSize;
      double Omega0;
      double OmegaLambda;
      double HubbleParam;
      char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];  /* fills to 256 Bytes */
    };


    struct io_header_3 //header for gadget3
    {
      int npart[6];
      /*!< number of particles of each type in this file */
      double mass[6];
      /*!< mass of particles of each type. If 0, then the masses are explicitly
                     stored in the mass-block of the snapshot file, otherwise they are omitted */
      double time;
      /*!< time of snapshot file */
      double redshift;
      /*!< redshift of snapshot file */
      int flag_sfr;
      /*!< flags whether the simulation was including star formation */
      int flag_feedback;
      /*!< flags whether feedback was included (obsolete) */
      unsigned int nPartTotal[6];
      /*!< total number of particles of each type in this snapshot. This can be
                       different from npart if one is dealing with a multi-file snapshot. */
      int flag_cooling;
      /*!< flags whether cooling was included  */
      int num_files;
      /*!< number of files in multi-file snapshot */
      double BoxSize;
      /*!< box-size of simulation in case periodic boundaries were used */
      double Omega0;
      /*!< matter density in units of critical density */
      double OmegaLambda;
      /*!< cosmological constant parameter */
      double HubbleParam;
      /*!< Hubble parameter in units of 100 km/sec/Mpc */
      int flag_stellarage;
      /*!< flags whether the file contains formation times of star particles */
      int flag_metals;
      /*!< flags whether the file contains metallicity values for gas and star
                     particles */
      unsigned int nPartTotalHighWord[6];
      /*!< High word of the total number of particles of each type (see header2)*/
      int flag_entropy_instead_u;
      /*!< flags that IC-file contains entropy instead of u */
      int flag_doubleprecision;
      /*!< flags that snapshot contains double-precision instead of single precision */

      int flag_ic_info;
      /*!< flag to inform whether IC files are generated with ordinary Zeldovich approximation,
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
      float lpt_scalingfactor;
      /*!< scaling factor for 2lpt initial conditions */

      char fill[48];    /*!< fills to 256 Bytes */

    };

/*!< holds header for snapshot files */



    template<typename FloatType>
    io_header_2 CreateGadget2Header(FloatType *masses, long *npart, double Boxlength,
                                    const cosmology::CosmologicalParameters<FloatType> &cosmology) {
      //types 2,3,4 are currently unused (only one zoom level)
      io_header_2 header2;
      header2.npart[0] = npart[0]; //gas
      header2.npart[1] = npart[1]; //higres
      header2.npart[2] = npart[2];
      header2.npart[3] = npart[3];
      header2.npart[4] = npart[4];
      header2.npart[5] = npart[5]; //lowres
      header2.mass[0] = masses[0];
      header2.mass[1] = masses[1];
      header2.mass[2] = masses[2];
      header2.mass[3] = masses[3];
      header2.mass[4] = masses[4];
      header2.mass[5] = masses[5];
      header2.time = cosmology.scalefactor;
      header2.redshift = cosmology.redshift;
      header2.flag_sfr = 0;
      header2.flag_feedback = 0;
      header2.nPartTotal[0] = npart[0];
      header2.nPartTotal[1] = npart[1];
      header2.nPartTotal[2] = npart[2];
      header2.nPartTotal[3] = npart[3];
      header2.nPartTotal[4] = npart[4];
      header2.nPartTotal[5] = npart[5];
      header2.flag_cooling = 0;
      header2.num_files = 1;
      header2.BoxSize = Boxlength;
      header2.Omega0 = cosmology.OmegaM0;
      header2.OmegaLambda = cosmology.OmegaLambda0;
      header2.HubbleParam = cosmology.hubble;

      if (npart[0] > 0) { //options for baryons
        header2.flag_sfr = 1;
        header2.flag_feedback = 1;
        header2.flag_cooling = 1;
      }


      return header2;
    }

    template<typename FloatType>
    io_header_3 CreateGadget3Header(FloatType *masses, long *npart, double Boxlength,
                                    const cosmology::CosmologicalParameters<FloatType> &cosmology) {

      //types 2,3,4 are currently unused (only one zoom level)

      io_header_3 header3;
      header3.npart[0] = (unsigned int) (npart[0]); //gas
      header3.npart[1] = (unsigned int) (npart[1]); //highres
      header3.npart[2] = (unsigned int) (npart[2]);
      header3.npart[3] = (unsigned int) (npart[3]);
      header3.npart[4] = (unsigned int) (npart[4]);
      header3.npart[5] = (unsigned int) (npart[5]); //lowres
      header3.mass[0] = (double) (masses[0]);
      header3.mass[1] = (double) (masses[1]);
      header3.mass[2] = (double) (masses[2]);
      header3.mass[3] = (double) (masses[3]);
      header3.mass[4] = (double) (masses[4]);
      header3.mass[5] = (double) (masses[5]);
      header3.time = cosmology.scalefactor;
      header3.redshift = cosmology.redshift;
      header3.flag_sfr = 0;
      header3.flag_feedback = 0;
      header3.nPartTotal[0] = (unsigned int) (npart[0]);
      header3.nPartTotal[1] = (unsigned int) (npart[1]);
      header3.nPartTotal[2] = (unsigned int) (npart[2]);
      header3.nPartTotal[3] = (unsigned int) (npart[3]);
      header3.nPartTotal[4] = (unsigned int) (npart[4]);
      header3.nPartTotal[5] = (unsigned int) (npart[5]);
      header3.flag_cooling = 0;
      header3.num_files = 1;
      header3.BoxSize = Boxlength;
      header3.Omega0 = cosmology.OmegaM0;
      header3.OmegaLambda = cosmology.OmegaLambda0;
      header3.HubbleParam = cosmology.hubble;
      header3.flag_stellarage = 0;  /*!< flags whether the file contains formation times of star particles */
      header3.flag_metals = 0;    /*!< flags whether the file contains metallicity values for gas and star  particles */
      header3.nPartTotalHighWord[0] = (unsigned int) (npart[0] >> 32);
      header3.nPartTotalHighWord[1] = (unsigned int) (npart[1] >> 32); //copied from Gadget3
      header3.nPartTotalHighWord[2] = (unsigned int) (npart[2] >> 32);
      header3.nPartTotalHighWord[3] = (unsigned int) (npart[3] >> 32);
      header3.nPartTotalHighWord[4] = (unsigned int) (npart[4] >> 32);
      header3.nPartTotalHighWord[5] = (unsigned int) (npart[5] >> 32);
      header3.flag_entropy_instead_u = 0; /*!< flags that IC-file contains entropy instead of u */
      header3.flag_doubleprecision = floatinfo<FloatType>::doubleprecision;
      header3.flag_ic_info = 1;
      header3.lpt_scalingfactor = 0.; /*!dummy value since we never use ic_info!=1 */

      if (npart[0] > 0) { //options for baryons & special behavior
        header3.flag_sfr = 1;
        header3.flag_feedback = 1;
        header3.flag_cooling = 1;
        //none of the following is currently supported in this code:
        //header3.flag_stellarage=1;  /*!< flags whether the file contains formation times of star particles */
        //header3.flag_metals=1;    /*!< flags whether the file contains metallicity values for gas and star  particles */
        //header3.flag_entropy_instead_u=1; /*!< flags that IC-file contains entropy instead of u */
      }

      return header3;
    }


#ifdef HAVE_HDF5
    template<typename FloatType>
void SaveHDF(const char* Filename, int n, io_header_3 header, FloatType *p_Data1,  FloatType *p_Data2, FloatType *p_Data3, FloatType *v_Data1,  FloatType *v_Data2, FloatType *v_Data3, const char *name1, const char *name2, const char *name3) //saves 3 arrays as datasets "name1", "name2" and "name3" in one HDF5 file
{
  hid_t       file_id, gidh, gid, dset_id,ndset_id,iddset_id, hdset_id,mt_id;   /* file and dataset identifiers */
  hid_t       nfilespace,dataspace_id, memspace, hboxfs,mt,idfilespace,IDmemspace,iddataspace_id;/* file and memory dataspace identifiers */
  hsize_t ncount[2],ncount2[2], noffset[2], h[1], stride[2],  idcount[2],idoffset[2],idstride[2],mta[1],idblock[2];
  h[0]=1;
  int m_nGrid[3]={n,n,n};
  mta[0]=6;

  file_id = H5Fcreate( Filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  gidh=H5Gcreate(file_id, "Header", 100 ); //100 gives max number of arguments (?)

   FloatType *box=new FloatType[1];
   FloatType *mtt=new FloatType[6];
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

   //hdset_id = H5Acreate( gidh, "NumPart_ThisFile", H5T_NATIVE_UINT, mt, H5P_DEFAULT );
   //H5Awrite( hdset_id, H5T_NATIVE_UINT, boxu);
   //H5Aclose(hdset_id);

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

   //boxi[0]=header.num_files;
   //hdset_id = H5Acreate( gidh, "NumFilesPerSnapshot", H5T_NATIVE_INT, hboxfs, H5P_DEFAULT );
   //H5Awrite( hdset_id, H5T_NATIVE_INT, boxi);
   //H5Aclose(hdset_id);


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

template<typename FloatType>
int save_phases(complex<FloatType> *phk, FloatType* ph, complex<FloatType> *delta, int n, const char *name){

    FloatType *helper=(FloatType*)calloc(n*n*n, sizeof(FloatType));
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

    template<typename FloatType>
    void SaveGadget2(const char *filename, long nPart, io_header_2 header1, FloatType *Pos1, FloatType *Vel1, FloatType *Pos2,
                     FloatType *Vel2, FloatType *Pos3, FloatType *Vel3, FloatType *Mass = NULL) {
      FILE *fd = fopen(filename, "w");
      if (!fd) throw std::runtime_error("Unable to open file for writing");
      FloatType *Pos = (FloatType *) calloc(3, sizeof(FloatType));
      int dummy;
      long i;
      //header block
      dummy = sizeof(header1);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      my_fwrite(&header1, sizeof(header1), 1, fd);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      //position block
      dummy = sizeof(FloatType) * (long) (nPart) *
              3; //this will be 0 or some strange number for n>563; BUT: gagdget does not actually use this value; it gets the number of particles from the header
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      for (i = 0; i < nPart; i++) {
        Pos[0] = Pos1[i];
        Pos[1] = Pos2[i];
        Pos[2] = Pos3[i];
        my_fwrite(Pos, sizeof(FloatType), 3, fd);
      }
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      //velocity block
      //dummy=sizeof(FloatType)*3*nPart; //this will be 0 or some strange number for n>563; BUT: gagdget does not actually use this value; it gets the number of particles from the header
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      for (i = 0; i < nPart; i++) {
        Pos[0] = Vel1[i];
        Pos[1] = Vel2[i];
        Pos[2] = Vel3[i];
        my_fwrite(Pos, sizeof(FloatType), 3, fd);
      }

      free(Pos);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      dummy = sizeof(long) *
              nPart; //here: gadget just checks if the IDs are ints or long longs; still the number of particles is read from the header file
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      for (i = 0; i < nPart; i++) {
        my_fwrite(&i, sizeof(long), 1, fd);

      }
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      if (Mass != NULL) {
        dummy = sizeof(FloatType) * nPart;
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
        my_fwrite(Mass, sizeof(FloatType), nPart, fd);
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
      }


      fclose(fd);

    }

    template<typename GridDataType, typename FloatType>
    void save(const std::string &name, double Boxlength,
              particle::mapper::ParticleMapper<GridDataType> & mapper,
              particle::AbstractMultiLevelParticleGenerator<GridDataType> & generator,
              const cosmology::CosmologicalParameters<FloatType> &cosmology, int gadgetformat) {

      std::cerr << "Hello from Gadget output!" << std::endl;

      FloatType tot_mass = 0.0;
      FloatType min_mass, max_mass, gas_mass;
      size_t nlow, ngas, nhigh;

      // input units [Gadget default]::
      // pmass in 1e10 h^-1 Msol
      // pos in Mpc h^-1
      // vel in km s^-1 a^1/2

      getParticleInfo(generator, mapper, min_mass, max_mass, tot_mass, gas_mass, ngas, nlow, nhigh);

      cerr << "min and max particle mass : " << min_mass << " " << max_mass << endl;
      cerr << "particle numbers (gas, high, low) : " << ngas << " " << nhigh << " " << nlow << endl;

      FloatType *masses = (FloatType *) calloc(6, sizeof(FloatType));
      long *npart = (long *) calloc(6, sizeof(long));

      //extend these arrays when additional particle types are added:
      masses[0] = gas_mass;
      masses[1] = min_mass;
      masses[5] = max_mass;
      npart[0] = ngas;
      npart[1] = nhigh;
      npart[5] = nlow;

      std::stringstream filename;
      filename << name << gadgetformat;
      FILE *fd;
      fd = fopen(filename.str().c_str(), "w");
      if (!fd) throw std::runtime_error("Unable to open file for writing");


      size_t n = npart[0] + npart[1] + npart[2] + npart[3] + npart[4] + npart[5];
      cerr << "total n = " << n << endl;
      int dummy;

      if (gadgetformat == 3) {
        io_header_3 header1 = CreateGadget3Header(masses, npart, Boxlength, cosmology);
        //header block
        dummy = sizeof(header1);
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
        my_fwrite(&header1, sizeof(header1), 1, fd);
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
      }
      else if (gadgetformat == 2) {
        io_header_2 header1 = CreateGadget2Header(masses, npart, Boxlength, cosmology);
        dummy = sizeof(header1);
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
        my_fwrite(&header1, sizeof(header1), 1, fd);
        my_fwrite(&dummy, sizeof(dummy), 1, fd);
      }
      else {
        cerr << "Wrong format for Gadget output!" << endl;
        exit(1);
      }

      float output_cache;

      dummy = sizeof(output_cache) * (long) (n) * 3;

      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      // DM particles positions

      for (auto i = mapper.beginDm(generator); i != mapper.endDm(generator); ++i) {
        auto particle = i.getParticle();

        output_cache = particle.pos.x;
        my_fwrite(&output_cache, sizeof(output_cache), 1, fd);
        output_cache = particle.pos.y;
        my_fwrite(&output_cache, sizeof(output_cache), 1, fd);
        output_cache = particle.pos.z;
        my_fwrite(&output_cache, sizeof(output_cache), 1, fd);

      }
      my_fwrite(&dummy, sizeof(dummy), 1, fd); //end of position block

      //cerr<< "counter after DM: "<< counter << endl;

      //gas particles velocities
      my_fwrite(&dummy, sizeof(dummy), 1, fd); //beginning of velocity block

      for (auto i = mapper.beginDm(generator); i != mapper.endDm(generator); ++i) {
        auto particle = i.getParticle();

        output_cache = particle.vel.x;
        my_fwrite(&output_cache, sizeof(output_cache), 1, fd);
        output_cache = particle.vel.y;
        my_fwrite(&output_cache, sizeof(output_cache), 1, fd);
        output_cache = particle.vel.z;
        my_fwrite(&output_cache, sizeof(output_cache), 1, fd);

      }
      my_fwrite(&dummy, sizeof(dummy), 1, fd); //end of velocity block

      //cerr<< "counter after DM vel: "<< counter << endl;

      //particle IDs (one for each gas, high res and low res particle)
      dummy = sizeof(long) *
              (n); //here: gadget just checks if the IDs are ints or long longs; still the number of particles is read from the header file
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
      for (size_t i = 0; i < n; i++) {
        my_fwrite(&i, sizeof(long), 1, fd);

        //counter+=1;
      }
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      //cerr<< "counter after IDs: "<< counter << endl;

      //IFF we want to save individual particle masses, they would go here, before the gas particle energies


      //gas particles energies: TODO TEST the prop. constant (especially unitv), (include YHe_ and gamma_ in the paramfile for baryons?)
      // if (npart[0] >0) {
      //     dummy=sizeof(FloatType)*(long)(n); //energy block
      //     my_fwrite(&dummy, sizeof(dummy), 1, fd);

      //     const FloatType YHe=0.248; //helium fraction, using default from MUSIC
      //     const FloatType gamma=5.0/3.0; //adiabatic gas index, using default from MUSIC

      //     const FloatType npol  = (fabs(1.0-gamma)>1e-7)? 1.0/(gamma-1.) : 1.0;
      //     const FloatType unitv = 1e5; //this is probably a unit transformation
      //     const FloatType h2    = hubble*hubble*0.0001;
      //     const FloatType adec  = 1.0/(160.*pow(Ob0*h2/0.022,2.0/5.0));
      //     const FloatType Tcmb0 = 2.726;
      //     const FloatType Tini  = a<adec? Tcmb0/a : Tcmb0/a/a*adec;
      //     const FloatType mu    = (Tini>1.e4) ? 4.0/(8.-5.*YHe) : 4.0/(1.+3.*(1.-YHe));
      //     FloatType ceint = 1.3806e-16/1.6726e-24 * Tini * npol / mu / unitv / unitv;

      //     cerr << "gas internal energy " << ceint << endl;

      //     for(auto i=pMapper->beginGas(); i!=pMapper->endGas(); ++i)
      //         my_fwrite(&ceint,sizeof(FloatType),1,fd);

      //     my_fwrite(&dummy, sizeof(dummy), 1, fd); //end of energy block
      // }



      fclose(fd);
      // free(Pos);

    }



  }
}


#endif //IC_GADGET_HPP
