"""
Script for converting .gadget3 binary ICs to the hdf5 format, for use with the EAGLE code.

Adapted from an IDL script by Rob Crain, by Jonathan Davies, 2021.

Usage: ics2hdf.py ics_file
* Creates an hdf5 file in the same directory as the original ICs, with the same name.
"""

import numpy as np
import os
import argparse
from scipy.io import FortranFile
import h5py as h5

def ics2hdf(ics_file):

    # Do quick check to make sure input is valid
    assert ics_file.name.split('.')[-1] == 'gadget3','Input ICs not valid. Please input GADGET3 binary ICs (.gadget3)'

    out_ics = ics_file.name[:-7] + 'hdf5'

    print(f"Reading binary file {ics_file.name}")

    # Open the FORTRAN unformatted binary ICs
    f = FortranFile(ics_file, 'r')

    # Read the header with the proper types
    header = f.read_record('(6,)i4','(6,)f8','f8','f8','i4','i4','(6,)i4','i4','i4','f8','f8','f8','f8','i4','i4','(6,)i4','i4','i4','(56,)b')

    hdict = {}
    hdict['NumPart_ThisFile'] = header[0]
    hdict['MassTable'] = header[1]
    hdict['ExpansionFactor'] = header[2][0]
    hdict['Time'] = header[2][0]
    hdict['Redshift'] = header[3][0]
    hdict['Flag_Sfr'] = header[4][0]
    hdict['Flag_Feedback'] = header[5][0]
    hdict['NumPart_Total'] = header[6]
    hdict['Flag_Cooling'] = header[7][0]
    hdict['NumFilesPerSnapshot'] = header[8][0]
    hdict['BoxSize'] = header[9][0]
    hdict['Omega0'] = header[10][0]
    hdict['OmegaLambda'] = header[11][0]
    hdict['HubbleParam'] = header[12][0]
    hdict['Flag_StellarAge'] = header[13][0]
    hdict['Flag_Metals'] = header[14][0]
    hdict['NumPart_Total_HighWord'] = header[15]
    hdict['Flag_DoublePrecision'] = header[17][0]

    # Correction to header (from Rob Crain's original IDL script):
    hdict['Flag_DoublePrecision'] = np.int32(0)

    # Pointer to this dict element for later convenience
    npart = hdict['NumPart_ThisFile']

    # Get the particle array lengths
    tot_npart = npart[1] + npart[2]

    # Load the particles
    pos = f.read_record('(%d,3)f4'%(tot_npart))
    vel = f.read_record('(%d,3)f4'%(tot_npart))
    ids = f.read_record('(%d,)i8'%(tot_npart))


    # Write out the ICs in hdf5 format

    print('Writing hdf5 file to '+out_ics)

    out = h5.File(out_ics,'w')

    print('Writing header...')

    out_header = out.create_group('Header')

    for key in hdict.keys():
        out_header.attrs.create(key,data=hdict[key])

    offset = 0
    for ptype in range(len(npart)):

        num = npart[ptype]
        if num == 0:
            continue

        print('Writing PartType'+str(ptype)+'...')

        out_pgroup = out.create_group('PartType%i'%(ptype))

        out_pgroup.create_dataset('Coordinates',data=pos[offset:offset+num,:])
        out_pgroup.create_dataset('Velocity',data=vel[offset:offset+num,:])
        out_pgroup.create_dataset('ParticleIDs',data=ids[offset:offset+num])

        offset += num

    print('Done!')



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'input_ICs',
        type=argparse.FileType(mode="br"),
        help="The input .gadget3 binary ICs file"
    )
    args = parser.parse_args()
    ics2hdf(args.input_ICs)


if __name__ == "__main__":
    main()
