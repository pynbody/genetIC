#!/usr/bin/env python

"""
Script to compare the output from two genetIC runs

Usage: compare.py path_to_output/
 * compares the particle output (path_to_output/*.gadget or path_to_output/*.tipsy) with path_to_output/reference_output
 * compares the grid output (path_to_output/grid-?.npy) with path_to_output/reference_grid
 * compares the power spectrum output (path_to_output/*.ps) with path_to_output/reference_ps/*.ps
"""


from __future__ import print_function
import numpy as np
import numpy.testing as npt
import os
import pynbody
import sys
import glob
import os.path
import warnings
import re
import platform

def infer_compare_decimal(sim):
    if sim['vel'].dtype==np.float64:
        return 5
    else:
        return 4

def compare(f1,f2) :
    assert f1.families() == f2.families()
    for fam in f1.families():
        assert len(f1[fam])==len(f2[fam])
    npt.assert_almost_equal(f1['mass'],f2['mass'],decimal=6)
    if 'eps' in f1.loadable_keys():
        npt.assert_almost_equal(f1['eps'],f2['eps'],decimal=6)

    compare_decimal = infer_compare_decimal(f1)

    try:
        npt.assert_almost_equal(f1['vel'],f2['vel'],decimal=compare_decimal)
    except:
        post_compare_diagnostic_plot(f1, f2)
        raise
    npt.assert_almost_equal(f1['pos'],f2['pos'],decimal=compare_decimal)
    if 'iord' in f1.loadable_keys():
        assert (f1['iord']==f2['iord']).all()
    print("Particle output matches")
    if 'overdensity' in f1.loadable_keys():
        npt.assert_almost_equal(f1['overdensity'],f2['overdensity'],decimal=5)
        print("Overdensity array output matches")

    if 'u' in f2.gas.loadable_keys():
        # NB intentional to check for u in f2's loadable_keys
        # older versions of genetIC did not write 'u' with gas output
        if not isinstance(f1, pynbody.snapshot.gadgethdf.GadgetHDFSnap):
            assert 'u' in f1.gas.loadable_keys()
            npt.assert_almost_equal(f1.gas['u'], f2.gas['u'], decimal=5)

def post_compare_diagnostic_plot(f1,f2):
    import pylab as p
    p.figure(figsize=(12,7))

    p.subplot(211)
    p.plot(f1['vx'], f2['vx']-f1['vx'], 'k,')
    p.xlabel("$v_x$")
    p.ylabel(r"$\Delta v_x$")

    p.subplot(212)
    p.plot(f1['x'], f2['x']-f1['x'], 'k,')
    p.xlabel("$x$")
    p.ylabel(r"$\Delta x$")

    p.savefig("diagnostic.png", bbox_inches='tight')
    p.close()
    _try_opening('diagnostic.png')

def _try_opening(filename):
    """Try opening the specified file in a gui"""
    if "Darwin" in platform.platform():
        os.system("open %s"%filename)
    elif "Linux" in platform.platform():
        os.system("xdg-open %s"%filename)

def compare_grids(ref, test):
    list_of_grids = [os.path.basename(x) for x in glob.glob(ref+"grid-?.npy")]
    list_of_grids.sort()
    assert (len(list_of_grids) != 0), "Could not find reference grids in the reference_grid folder"
    for grid in list_of_grids:
        grid_ref = np.load(ref+grid)
        grid_test = np.load(test+grid)
        try:
            npt.assert_almost_equal(grid_ref, grid_test, decimal=4)
        except:
            post_compare_grid_diagnostic_plot(ref, test)
            raise
    print("Grid output matches")

def post_compare_grid_diagnostic_plot(ref, test):
    import plotslice as ps
    import pylab as p
    p.figure(figsize=(12,7))
    p.set_cmap("RdBu_r")
    p.subplot(131)
    ps.plotslice(ref)
    p.title("Reference")
    p.subplot(132)
    ps.plotslice(test)
    p.title("Test result")
    p.subplot(133)
    p.title("Diff")
    ps.plotslice(test, diff_prefix=ref)
    p.savefig("diagnostic_grids.pdf", bbox_inches='tight')
    p.close()
    _try_opening('diagnostic_grids.pdf')

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

def _strip_time(line):
    return re.sub("^ [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}   ", "", line)

def compare_outputlogfile(reference_file, test_file):
    with open(test_file) as f:
        with open(reference_file) as ref:
            lines_to_match = ref.read().splitlines()
            lines_in_output = f.read().splitlines()
            lines_in_output = [_strip_time(t) for t in lines_in_output]
            for s in lines_to_match:
                assert(s in lines_in_output), "Line %s from reference.txt is not present in the logged output" % s
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
    gadget_multifiles = glob.glob(dirname + "/*.gadget?.0")
    if len(gadget_multifiles) == 1:
        gadget_multifiles = [f[:-2] for f in gadget_multifiles]

    hdf_multifiles = glob.glob(dirname + "/*.0.hdf5")
    if len(hdf_multifiles) == 1:
        hdf_multifiles = [f[:-7] for f in hdf_multifiles]

    hdf_files = glob.glob(dirname + "/*.hdf5")
    hdf_files = [f for f in hdf_files if not any([f.startswith(h) for h in hdf_multifiles])]

    return (glob.glob(dirname + "/*.tipsy") + glob.glob(dirname + "/*.gadget?") + gadget_multifiles +
            hdf_files + hdf_multifiles)


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
    if len(output_file)>1:
        raise RuntimeError("More than one particle file to test against. Don't know which to use")
    if len(output_file)>0 and os.path.exists(sys.argv[1]+"/reference_output"):
        compare(pynbody.load(output_file[0]),pynbody.load(sys.argv[1]+"/reference_output"))

if __name__=="__main__":
    warnings.simplefilter("ignore")
    if len(sys.argv)==2:
        default_comparisons()
    else:
        compare(pynbody.load(sys.argv[1]), pynbody.load(sys.argv[2]))
