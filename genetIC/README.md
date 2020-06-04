Very brief installation instructions
------------------------------------

You need a modern C++ compiler, `fftw3` and `GSL`. For mac, you can get these 
from the macports project. I recommend installing packages `gcc9`, `fftw-3` 
and `gsl`. With these installed, you should be able to simply `make` the project. 

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
