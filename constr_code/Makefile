HOST	= $(shell if [ `which scutil` ]; then scutil --get ComputerName | cut -c1-5; else hostname | cut -c1-5 ; fi )

CC    	= g++
C	= gcc
CFLAGS	= -Wall -g -fopenmp -std=c++11   -fdiagnostics-color=auto
CPATH = /opt/local/include/
LPATH   =/opt/local/lib


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


ifeq ($(HOST), Piann)
	CPATH = /sw/include/
	LPATH = /sw/lib/
	C = gcc-4
	CC = g++-4
endif

ifeq ($(HOST), Rhodo)
	C = clang
	CC = clang -lstdc++
endif

ifeq ($(HOST), dirac)
	CC 	= icc
	FFTW    = -I$(GSL_INC_DIR)
	FFTWLIB = -L$(FFTWDIR) -I$(FFTWINCLUDE)  -lfftw -lrfftw
endif

ifeq ($(HOST), chewb)
	CC      = /usr/local/bin/g++
	FFTW    = -DFFTW_TYPE_PREFIX=1
            FFTWLIB = -ldfftw -ldrfftw
endif


#HDF5_FLAGS = -lhdf5 -D H5_USE_16_API -DHAVE_HDF5

all: IC


IC: src/main.cpp src/HDF_IO.hh src/ic.hpp src/constrainer.hpp src/cosmo.hpp src/fft.hpp src/float_types.hpp src/grid.hpp src/io.hpp src/parser.hpp src/sparse.hpp Makefile
	$(CC) $(CFLAGS) -o IC -I$(CPATH) $(FFTW) src/main.cpp -L$(LPATH) -lm $(FFTWLIB) -lgsl -lgslcblas $(HDF5_FLAGS) -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
clean:
	rm -rf IC
