CC    	= g++
C	= gcc
CFLAGS	= -Wall -O2
CPATH = /usr/include
LPATH   =/usr/local/lib
HDF5_FLAGS = -lhdf5 -D H5_USE_16_API -DHAVE_HDF5

all: IC_deriv

IC_deriv: ICgauss_deriv_clean.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_deriv -I$(CPATH) ICgauss_deriv_clean.cc -L$(LPATH) -lm -ldrfftw -ldfftw -lgsl -lgslcblas $(HDF5_FLAGS) -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION 
clean: 
	rm -rf IC_deriv
