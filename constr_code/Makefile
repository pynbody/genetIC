HOST	= $(shell if [ `which scutil` ]; then scutil --get ComputerName | cut -c1-5; else hostname | cut -c1-5 ; fi )
NUMPY_INCLUDE = -Ilibnpy-0.5/include/ -Llibnpy-0.5/lib/ -lnpy

CC    	= g++
C	= gcc
CFLAGS	= -Wall -g -fopenmp -std=c++11   -fdiagnostics-color=auto
CPATH = /opt/local/include/
LPATH   =/opt/local/lib

ifeq ($(HOST), Piann)
	CPATH = /sw/include/
	LPATH = /sw/lib/
	C = gcc-4
	CC = g++-4
endif

ifeq ($(HOST), dirac)
	CC 	= icc
	FFTW    = -I$(GSL_INC_DIR)
	FFTWLIB = -L$(FFTWDIR) -I$(FFTWINCLUDE)  -lfftw -lrfftw
else
	ifeq ($(HOST), chewb)
		CC      = /usr/local/bin/g++
		FFTW    = -DFFTW_TYPE_PREFIX=1
                FFTWLIB = -ldfftw -ldrfftw
	else
		FFTW    = -DFFTW_TYPE_PREFIX=1
		FFTWLIB = -ldfftw -ldrfftw
	endif
endif
#HDF5_FLAGS = -lhdf5 -D H5_USE_16_API -DHAVE_HDF5

all: IC


IC: src/main.cpp src/HDF_IO.hh src/ic.hpp src/constrainer.hpp src/cosmo.hpp src/fft.hpp src/float_types.hpp src/grid.hpp src/io.hpp src/parser.hpp src/sparse.hpp Makefile
	$(CC) $(CFLAGS) $(NUMPY_INCLUDE) -o IC -I$(CPATH) $(FFTW) src/main.cpp -L$(LPATH) -lm $(FFTWLIB) -lgsl -lgslcblas $(HDF5_FLAGS) -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
clean:
	rm -rf IC
