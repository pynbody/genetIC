# Reasonably sensible defaults
CC    	= g++
C		= gcc
CFLAGS	= -Wall -g -O0 -fopenmp -std=c++14   -fdiagnostics-color=auto
CPATH 	= /opt/local/include/
LPATH   = /opt/local/lib


# Use -DFFTW3 to use FFTW3 (instead of FFTW 2).
#
#
# Additionally use -DFFTW3_THREADS=n to use threads. This will use n threads
# UNLESS you are compiling with openmp, in which case you do not need to set n
# as it will automatically be set to omp_get_num_threads()
#
# Note you have to link to fftw3 and fftw3_threads (or fftw3_omp)
# if you are using threads

FFTW = -DFFTW3 -DFFTW_THREADS
FFTWLIB = -lfftw3 -lfftw3_threads

# Example for FFTW2:

#FFTW    = -DFFTW_TYPE_PREFIX=1
#FFTWLIB = -ldfftw -ldrfftw


# If you want to use HDF5, something like the following is necessary.
#
#HDF5_FLAGS = -lhdf5 -D H5_USE_16_API -DHAVE_HDF5
#
# WARNING - HDF5 has not been tested recently and in particular won't work
# with zoom simulations as the code does not support varying mass particles




# Some customizations for known systems

HOST	= $(shell if [ `which scutil 2>/dev/null` ]; then scutil --get ComputerName | cut -c1-5; else hostname | cut -c1-5 ; fi )
GIT_VERSION = $(shell git rev-parse --short HEAD)
GIT_MODIFIED = $(shell git ls-files -m | tr "\n" " ")
GIT_VARIABLES = -DGIT_VERSION='"$(GIT_VERSION)"' -DGIT_MODIFIED='"$(GIT_MODIFIED)"'


ifeq ($(HOST), bridg)
	CC = /nasa/pkgsrc/2015Q1/gcc49/bin/gcc
	CPATH = /nasa/gsl/1.14/include/:/nasa/intel/Compiler/2015.0.090/composer_xe_2015.0.090/mkl/include/fftw/:/nasa/intel/Compiler/2015.0.090/composer_xe_2015.0.090/mkl/include/
	LPATH =  /nasa/gsl/1.14/lib/:/nasa/intel/Compiler/2015.0.090/composer_xe_2015.0.090/mkl/lib/intel64/
	FFTWLIB = -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -lstdc++ -lgfortran -lfftw3xc_intel
	FFTW = -DFFTW3 # not working with threads for some reason :-(
endif

ifeq ($(HOST), Piann)
	CPATH = /sw/include/
	LPATH = /sw/lib/
	C = gcc-4
	CC = g++-4
	FFTW    =
	FFTWLIB =  -lfftw -lrfftw
endif

ifeq ($(HOST), Rhodo)
	C = /usr/bin/g++
	CC = /usr/bin/g++ -lstdc++ -O3 -g -L/usr/local/lib/
endif

ifeq ($(HOST), dirac)
	CC 	= g++ -O3 -g
	FFTW    = -I$(GSL_INC_DIR)
	FFTWLIB = -L$(FFTWDIR) -I$(FFTWINCLUDE)  -lfftw -lrfftw
endif

ifeq ($(HOST), chewb)
	CC      = /usr/local/bin/g++
	FFTW    = -DFFTW_TYPE_PREFIX=1
    FFTWLIB = -ldfftw -ldrfftw
endif



all: IC


IC: src/main.cpp src/filesystem.cpp src/HDF_IO.hh src/ic.hpp src/constraint.hpp src/multilevelfieldmanager.hpp src/cosmo.hpp src/fft.hpp src/float_types.hpp src/grid.hpp src/io.hpp src/parser.hpp src/sparse.hpp src/mapper.hpp src/filter.hpp Makefile
	$(CC) $(CFLAGS) -o IC $(GIT_VARIABLES) -I$(CPATH) $(FFTW) src/main.cpp src/filesystem.cpp -L$(LPATH) -lm $(FFTWLIB) -lgsl -lgslcblas $(HDF5_FLAGS) -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
clean:
	rm -rf IC
