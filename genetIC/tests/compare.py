#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import numpy.testing as npt
import pynbody
import sys
import glob
import os.path
import warnings

def compare(f1,f2) :
    npt.assert_almost_equal(f1['mass'],f2['mass'],decimal=6)
    if 'eps' in f1.loadable_keys():
        npt.assert_almost_equal(f1['eps'],f2['eps'],decimal=6)

    if f1['vel'].dtype==np.float64:
        compare_decimal = 5
    else:
        compare_decimal = 4

    npt.assert_almost_equal(f1['vel'],f2['vel'],decimal=compare_decimal)
    npt.assert_almost_equal(f1['pos'],f2['pos'],decimal=compare_decimal)
    print("Particle output matches")
    if 'overdensity' in f1.loadable_keys():
        npt.assert_almost_equal(f1['overdensity'],f2['overdensity'],decimal=5)
        print("Overdensity array output matches")

def compare_grids(ref, test):
    list_of_grids = [os.path.basename(x) for x in glob.glob(ref+"grid-?.npy")]
    list_of_grids.sort()
    assert (len(list_of_grids) != 0), "Could not find reference grids in the reference_grid folder"
    for grid in list_of_grids:
        grid_ref = np.load(ref+grid)
        grid_test = np.load(test+grid)
        npt.assert_almost_equal(grid_ref, grid_test, decimal=4)
    print("Grid output matches")

def compare_ps(ref, test):
    ref_vals = np.loadtxt(ref)
    test_vals = np.loadtxt(test)
    npt.assert_allclose(ref_vals, test_vals, rtol=1e-4)
    print("Power-spectrum output %s matches" % ref)


def check_comparison_is_possible(dirname):
    # A valid test must have either a tipsy output and its reference output or numpy grids and their references.

    output_file = particle_files_in_dir(dirname)
    assert(len(output_file)>=0)

    output_grids = [os.path.basename(x) for x in glob.glob(dirname+"grid-?.npy")]
    assert(len(output_grids)>=0)

    if len(output_file)==0:
    # If a tipsy test output is not generated, you must at least provide numpy grids
        if len(output_grids)==0:
            raise IOError("There are no particle files or numpy files to test against")
        elif not os.path.exists(dirname+"/reference_grid"):
            raise IOError("You must provide a reference_grid folder to test against if no particle output is generated")

    elif len(output_file)==1:
        if not os.path.isfile(sys.argv[1]+"/reference_output"):
            raise IOError("A particle output is present but there is no reference_output to test against.")

    else:
        raise IOError("There is more than one particle output file to test against.")


def particle_files_in_dir(dirname):
    return glob.glob(dirname + "/*.tipsy") + glob.glob(dirname + "/*.gadget?")


def default_comparisons():

    assert len(sys.argv)==2
    check_comparison_is_possible(sys.argv[1])

    if os.path.exists(sys.argv[1]+"/reference_grid"):
        compare_grids(sys.argv[1]+"/reference_grid/",sys.argv[1]+"/")

    powspecs = sorted(glob.glob(sys.argv[1]+"/*.ps"))
    powspecs_test = sorted(glob.glob(sys.argv[1]+"/reference_ps/*.ps"))
    for ps, ps_test in zip(powspecs, powspecs_test):
        compare_ps(ps,ps_test)

    output_file = particle_files_in_dir(sys.argv[1])
    if len(output_file)>0:
        compare(pynbody.load(output_file[0]),pynbody.load(sys.argv[1]+"/reference_output"))

if __name__=="__main__":
    warnings.simplefilter("ignore")
    if len(sys.argv)==2:
        default_comparisons()
    else:
        compare(pynbody.load(sys.argv[1]), pynbody.load(sys.argv[2]))
