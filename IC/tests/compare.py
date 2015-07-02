#!/usr/bin/env python

import numpy as np
import numpy.testing as npt
import pynbody
import sys

def compare(f1,f2) :
    npt.assert_almost_equal(f1['vel'],f2['vel'],decimal=3)
    npt.assert_almost_equal(f1['pos'],f2['pos'],decimal=4)
    npt.assert_almost_equal(f1['mass'],f2['mass'],decimal=6)
    #npt.assert_almost_equal(f1['eps'],f2['eps'])

if __name__=="__main__":
    assert len(sys.argv)==3
    compare(pynbody.load(sys.argv[1]),pynbody.load(sys.argv[2]))
