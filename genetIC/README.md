Very brief installation instructions
------------------------------------

The development team for `GenetIC` have tested the code on Linux and MacOS. 
Using any other operating system is not recommended, but may work (e.g. with WSL). 
You need a modern C++ compiler supporting OpenMP, plus the libraries `fftw3` and `GSL`. 
For MacOS, you can get these from the [macports](https://www.macports.org) project. 
I recommend installing packages `gcc9`, `fftw-3` and `gsl`. With these installed, 
you should be able to simply `make` the project. 

For more detailed installation instructions and explanation, see the PDF user
manual. This is available alongside the latest release at 
https://github.com/pynbody/genetIC/releases


Test suite
----------

The test suite requires python and pynbody to be installed. The easiest way
to accomplish this is to install `anaconda` python on your machine, then
type `pip install pynbody`.

Once compiled, enter the `tests` subfolder and type `./run_tests.sh`.
If all is OK, you will see `Tests seem OK`. You should also test the particle
mapping by typing `./run_mapper_tests.sh` which should also report that 
`Tests seem OK`.

For more information, see the PDF user manual at 
https://github.com/pynbody/genetIC/releases.
