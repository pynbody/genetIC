HOST	= $(shell hostname | cut -c1-5)

CC    	= g++
C	= gcc
CFLAGS	= -Wall -g -fopenmp -std=c++11   -fdiagnostics-color=auto
CPATH = /opt/local/include/
LPATH   =/opt/local/lib

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

all: IC_deriv

IC_deriv: main.cpp ICgauss_deriv_clean.cc Makefile HDF_IO.hh constrainer.h
	$(CC) $(CFLAGS) -o IC_deriv -I$(CPATH) $(FFTW) main.cpp -L$(LPATH) -lm $(FFTWLIB) -lgsl -lgslcblas $(HDF5_FLAGS) -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
clean:
	rm -rf IC_deriv
