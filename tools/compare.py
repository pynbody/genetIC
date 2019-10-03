#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import numpy.testing as npt
import pynbody
import sys
import glob
import os.path
import warnings
import re

def infer_compare_decimal(sim):
    if sim['vel'].dtype==np.float64:
        return 5
    else:
        return 4

def compare(f1,f2) :
    npt.assert_almost_equal(f1['mass'],f2['mass'],decimal=6)
    if 'eps' in f1.loadable_keys():
        npt.assert_almost_equal(f1['eps'],f2['eps'],decimal=6)

    compare_decimal = infer_compare_decimal(f1)

    npt.assert_almost_equal(f1['vel'],f2['vel'],decimal=compare_decimal)
    npt.assert_almost_equal(f1['pos'],f2['pos'],decimal=compare_decimal)
    if 'iord' in f1.loadable_keys():
        assert (f1['iord']==f2['iord']).all()
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

def compare_masks(ref, test):
    list_of_masks = [os.path.basename(x) for x in glob.glob(ref+"mask-?.npy")]
    list_of_masks.sort()
    assert (len(list_of_masks) != 0), "Could not find reference grids in the reference_grid folder"
    for mask in list_of_masks:
        mask_ref = np.load(ref+mask)
        mask_test = np.load(test+mask)
        npt.assert_almost_equal(mask_ref, mask_test, decimal=4)
    print("Mask output matches")

def compare_ps(ref, test):
    ref_vals = np.loadtxt(ref)
    test_vals = np.loadtxt(test)
    npt.assert_allclose(ref_vals, test_vals, rtol=1e-4)
    print("Power-spectrum output %s matches" % ref)

def compare_outputlogfile(reference_file, test_file):
    with open(test_file) as f:
        with open(reference_file) as ref:
            list = ref.read().splitlines()
            tests = f.read().splitlines()
            for s in list:
                assert(s in tests), "Line %s from reference.txt is not present in the logged output" % s
    print("Log output matches")

def get_grafic_files(path):
    return sorted(glob.glob(os.path.join(path,"*.grafic_*")))

def get_grafic_resolutions(files):
    res = [int(re.match(r".*\.grafic_([0-9]*)$", f).group(1)) for f in files]
    return res

def compare_grafic(reference_path, test_path):
    ref_files = get_grafic_files(reference_path)
    test_files = get_grafic_files(test_path)
    ref_res = get_grafic_resolutions(ref_files)
    test_res = get_grafic_resolutions(test_files)
    assert ref_res==test_res, "Grafic output files do not match in resolutions available"

    for fname1,fname2,res in zip(test_files, ref_files, ref_res):
        print("Test resolution",res)
        f1 = pynbody.load(fname1)
        f2 = pynbody.load(fname2)
        compare(f1,f2)
        assert (f1['iord']==f2['iord']).all()
        npt.assert_almost_equal(f1['deltab'],f2['deltab'],decimal=infer_compare_decimal(f1))

def check_comparison_is_possible(dirname):
    # A valid test must have either a tipsy/gadget output and its reference output or numpy grids and their references.

    if os.path.exists(dirname+"/reference.txt"):
        return # OK if we are just looking at the textual output

    output_file = particle_files_in_dir(dirname)
    assert(len(output_file)>=0)

    output_grids = [os.path.basename(x) for x in glob.glob(dirname+"grid-?.npy")]
    assert(len(output_grids)>=0)

    output_grafic = get_grafic_files(dirname)

    if len(output_file)==0 and len(output_grafic)==0:
        # If a particle test output is not generated, you must at least provide numpy grids
        if len(output_grids)==0:
            raise IOError("There are no particle files or numpy files to test against")
        elif not os.path.exists(dirname+"/reference_grid"):
            raise IOError("You must provide a reference_grid folder to test against if no particle output is generated")

    elif len(output_grafic)>0 and len(output_file)==0:
        pass # no specific tests implemented here
    elif len(output_file)==1 and len(output_grafic)==0:
        if not os.path.isfile(sys.argv[1]+"/reference_output"):
            raise IOError("A particle output is present but there is no reference_output to test against.")

    else:
        raise IOError("There is more than one particle output file to test against.")


def particle_files_in_dir(dirname):
    return glob.glob(dirname + "/*.tipsy") + glob.glob(dirname + "/*.gadget?")


def default_comparisons():

    assert len(sys.argv)==2
    check_comparison_is_possible(sys.argv[1])

    if os.path.exists(sys.argv[1]+"/reference.txt"):
        compare_outputlogfile(sys.argv[1]+"/reference.txt", sys.argv[1]+"/IC_output.txt")

    if os.path.exists(sys.argv[1]+"/reference_grid"):
        compare_grids(sys.argv[1]+"/reference_grid/",sys.argv[1]+"/")

    if os.path.exists(sys.argv[1]+"/reference_mask"):
        compare_masks(sys.argv[1]+"/reference_mask/",sys.argv[1]+"/")

    if os.path.exists(sys.argv[1]+"/reference_grafic/"):
        compare_grafic(sys.argv[1]+"/reference_grafic/", sys.argv[1])

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
