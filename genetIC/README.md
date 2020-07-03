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

The test suite requires Python 3 and `pynbody` to be installed. The easiest way
to accomplish this is to install [Anaconda Python 3](https://anaconda.org) on your
machine, then type `pip install pynbody`.

Once compiled, enter the `tests` subfolder and type `./run_tests.sh`.
If all is OK, you will see `Tests seem OK`. You should also test the particle
mapping by typing `./run_mapper_tests.sh` which should also report that 
`Tests seem OK`.

For more information, see the PDF user manual at 
https://github.com/pynbody/genetIC/releases.

Using Docker
------------

If you have [docker](https://docker.com) installed, you can conveniently use genetIC inside a container,
along with pynbody for running tests. 

To get the latest version from dockerhub, simply type
```
docker pull apontzen/genetic
```
*Alternatively*, to build your own docker image from this folder type

```
docker build -t apontzen/genetic .
```

Either of the two commands above creates an image and and are reasonably self-explanatory.

After the build or pull is complete, to run the example, type
```
cd ../example
docker run --rm -v `pwd`:/working_folder/ apontzen/genetic /working_folder/paramfile.txt  
```

To explain what is going on here:

 * `docker run` is the command to run a docker instance
 * `--rm` specifies the container should be removed on exit (don't worry, this doesn't delete the output data which is
 safely stored in the working directory)
 * ``-v `pwd`:/working_folder/`` mounts the example folder into the docker image filesystem, at
  `/working_folder/`. 
 * `apontzen/genetic` specifies the image you just created at the pull or build step 
 * `/working_folder/paramfile.txt` is the argument to `genetIC` that specifies the initial conditions parameter file.
 
To run the tests, use
```
cd ..
<<<<<<< HEAD
docker run --rm -v `pwd`:/genetic_repository/ -e IC=/genetIC/genetIC -w /genetic_repository/genetIC/tests --entrypoint bash apontzen/genetic run_tests.sh
docker run --rm -v `pwd`:/genetic_repository/ -e IC=/genetIC/genetIC -w /genetic_repository/genetIC/tests --entrypoint bash apontzen/genetic run_mapper_tests.sh
=======
docker run -v `pwd`:/genetic_repository/ -e IC=/genetIC/genetIC -w /genetic_repository/genetIC/tests genetic bash run_tests.sh
docker run -v `pwd`:/genetic_repository/ -e IC=/genetIC/genetIC -w /genetic_repository/genetIC/tests genetic bash run_mapper_tests.sh
>>>>>>> 7f29693d65e56db580292e6da5e15efd8d558241
```

