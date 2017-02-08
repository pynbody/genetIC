Pre-requisites:
---------------

You need `fftw3` and `GSL`. For mac, you can get these from the macports
project. I recommend installing packages `gcc5`, `fftw-3` and `gsl`.

You may need to edit the Makefile for your machine. Don't edit the global
settings in the makefile, create a new section like the existing ones for
`pfe`, `dirac` (which are typical Linux configurations) and `Snowd`, `Rhodo`
(which are typical Mac OS configurations).

Then to make, just type `make`.


Test suite
----------

The test suite requires python and pynbody to be installed. The easiest way
to accomplish this is to install `anaconda` python on your machine, then
type `pip install pynbody`.

Once compiled, enter the `tests` subfolder and type `./run_tests.sh`.
If all is OK, you will see `Tests seem OK`.
