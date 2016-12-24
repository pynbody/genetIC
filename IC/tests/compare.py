#!/usr/bin/env python

import numpy as np
import numpy.testing as npt
import pynbody
import sys
import glob
import os.path

def compare(f1,f2) :
    npt.assert_almost_equal(f1['vel'],f2['vel'],decimal=4)
    npt.assert_almost_equal(f1['pos'],f2['pos'],decimal=4)
    npt.assert_almost_equal(f1['mass'],f2['mass'],decimal=6)
    npt.assert_almost_equal(f1['eps'],f2['eps'],decimal=6)

def compare_grids(ref, test):
    list_of_grids = [os.path.basename(x) for x in glob.glob(ref+"grid-?.npy")]
    list_of_grids.sort()
    for grid in list_of_grids:
        grid_ref = np.load(ref+grid)
        grid_test = np.load(test+grid)
        npt.assert_almost_equal(grid_ref, grid_test, decimal=4)

if __name__=="__main__":
    assert len(sys.argv)==2
    if os.path.exists(sys.argv[1]+"/reference_grid"):
        compare_grids(sys.argv[1]+"/reference_grid/",sys.argv[1]+"/")
    output_file = glob.glob(sys.argv[1]+"/*.tipsy")
    assert len(output_file)==1, "Could not find a unique output file to test against"
    compare(pynbody.load(output_file[0]),pynbody.load(sys.argv[1]+"/reference_output"))
