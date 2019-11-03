genetIC python tools
--------------------

This folder contains python tools related to `genetIC`, chiefly for development purposes. These are compatible with python 3.5 and higher. Python 2.7 is _not_ supported.
You also need to install numpy, scipy, matplotlib and [pynbody](http://pynbody.github.io/pynbody/).

* `compare.py` is a module for comparing the output of a `genetIC` run with reference outputs, used by the test suite in `genetIC/tests`.
* `plotslice.py` is a module for plotting 2D and 1D slices of the output (either particle or grids).
* `check_zoom_accuracy.py` is a module for generating comparisons between zoom simulations and idealised counterparts (i.e. where the whole box is realised at the highest resolution). When run from the command line, it auto-generates a number of comparison tests and places the resulting figures in a `figures/`.

