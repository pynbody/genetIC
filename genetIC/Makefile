# Reasonably sensible defaults
CFLAGS	= -Wall -g -lpthread -O2 -fopenmp -std=c++14   -fdiagnostics-color=auto -I`pwd`

# Here are code options. Probably don't fiddle with these unless you know what you're doing.
#
# -DCUBIC_INTERPOLATION: Use cubic interpolation when moving information from a coarse to a fine grid. Slower but
#    much better.
# -DVELOCITY_MODIFICATION_GRADIENT_FOURIER_SPACE: As the name suggests, use "ik" as the gradient operator; otherwise
#    uses a fourth order stencil in real space. Enabling this is essential if you disable CUBIC_INTERPOLATION.
# -DFRACTIONAL_K_SPLIT: This defines the characteristic k0 of the filter used to split information between zoom levels,
#    defined by k0 = k_nyquist * FRACTIONAL_K_SPLIT, where k_nyquist is the Nyquist frequency of the lower-resolution grid.
#    Higher values correlate the zoom grids more precisely with the unzoomed base grids, but at the cost of errors
#    on small scales. The default value is an empirically-determined trade-off which works well.
# -DFRACTIONAL_FILTER_TEMPERATURE: This defines the characteristic temperature of the Fermi filter used to split
#    information between levels, as k0 * FRACTIONAL_FILTER_TEMPERATURE. Smaller values give sharper transitions in the
#    filtering, which leads to better accuracy in Fourier space but worse accuracy in real space. The default value is
#    an empirically-determined trade-off which works well.
# -DFILTER_ON_COARSE_GRID: When moving long wavelength information onto a zoom grid, filters on the base grid
#    and then interpolates into the zoom, rather than the default which is to interpolate and then filter.

CODEOPTIONS  = -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -DCUBIC_INTERPOLATION -DFRACTIONAL_K_SPLIT=0.3 -DFRACTIONAL_FILTER_TEMPERATURE=0.1 -DFILTER_ON_COARSE_GRID
CPATH 	= /opt/local/include/
LPATH   = /opt/local/lib

GSLFLAGS = -lgsl -lgslcblas

# Use -DFFTW3 to use FFTW3 (instead of FFTW 2).
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


# Provide information on the git repository, to be printed out on each run
GIT_MODIFIED = $(shell git ls-files -m 2>/dev/null | tr "\n" " " )
GIT_VERSION = $(shell git rev-parse --short HEAD 2>/dev/null)
GIT_VARIABLES = -DGIT_VERSION='"$(GIT_VERSION)"' -DGIT_MODIFIED='"$(GIT_MODIFIED)"'


# Some customizations for known systems

HOST	= $(shell if [ `which scutil 2>/dev/null` ]; then scutil --get ComputerName | cut -c1-5; else hostname | cut -c1-5 ; fi )
HOST3	= $(shell hostname | cut -c1-3)

ifeq ($(HOST3), pfe)
	CXX = /nasa/pkgsrc/2016Q2/gcc5/bin/g++
	CPATH = /u/apontzen/genetIC/genetIC:/nasa/gsl/1.14/include/:/nasa/intel/Compiler/2016.2.181/compilers_and_libraries_2016.2.181/linux/mkl/include/fftw
	LPATH =  /nasa/gsl/1.14/lib/ -L/nasa/intel/Compiler/2016.2.181/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64/
	CFLAGS = -Wall -g -O3 -fopenmp -std=c++14 -fdiagnostics-color=auto -I`pwd` -DOPENMP -DDOUBLEPRECISION
	FFTWLIB = -ldl -lpthread -lm -lstdc++ -lgfortran -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core
	FFTW = -DFFTW3 -DFFTW_THREADS
	GSLFLAGS = /nasa/gsl/1.14/lib/libgsl.a /nasa/gsl/1.14/lib/libgslcblas.a
endif

ifeq ($(HOST), Snowd)
	CPATH = /opt/local/include/
	LPATH = /opt/local/lib/
	CFLAGS = -O3 -fopenmp -std=c++14 -fdiagnostics-color=auto -I`pwd` -DOPENMP -DDOUBLEPRECISION
	CXX = /opt/local/bin/clang++-mp-6.0
endif

ifeq ($(HOST), clema)
        CPATH = /opt/local/include/
        LPATH = /opt/local/lib/
	      CFLAGS = -O2 -g -fopenmp -std=c++14 -fdiagnostics-color=auto -I`pwd` -DOPENMP -DDOUBLEPRECISION
        CXX = /opt/local/bin/clang++-mp-6.0
 endif

ifeq ($(HOST), jasmi)
	CXX = /opt/local/bin/g++-mp-9   -L/opt/local/lib/
	CFLAGS = -O2 -g -fopenmp -std=c++14 -I`pwd` -DOPENMP -DDOUBLEPRECISION
endif

ifeq ($(HOST), Rhodo)
	CXX = /usr/local/bin/g++ -lstdc++ -fopenmp -O3 -g -L/usr/local/lib/
endif

ifeq ($(HOST), dirac)
	CXX	= g++
	CFLAGS  += -m64 -I${FFTWINCLUDE} -I${GSL_INC_DIR}
        FFTW    = -DFFTW3
        FFTWLIB = -lgomp -lpthread -lm -ldl -lfftw3 -L${FFTWDIR} -L${GSL_LIB_DIR}
endif

ifeq ($(HOST), chewb)
	CXX      = /usr/local/bin/g++
	FFTW    = -DFFTW_TYPE_PREFIX=1
	FFTWLIB = -ldfftw -ldrfftw
endif

ifeq ($(HOST), imac7)
	CXX      = /opt/local/bin/g++-mp-6
	CFLAGS  += -Wextra
	FFTW    = -DFFTW3 -DFFTW_THREADS
	FFTWLIB = -lfftw3 -lfftw3_threads
endif

ifeq ($(HOST), marti)
	CXX      = /opt/local/bin/g++-mp-6
	CFLAGS  += -Wextra
	FFTW    = -DFFTW3 -DFFTW_THREADS
	FFTWLIB = -lfftw3 -lfftw3_threads
endif

ifeq ($(HOST3), hyp)
        CXX = g++
        CPATH = /share/apps/fftw-3.3.5/include/:/share/apps/gsl-2.4/include
        LPATH = /share/apps/gsl-2.4/lib
        LD_LIBRARY_PATH = /share/apps/gsl-2.4/lib:/share/apps/fftw-3.3.5/lib:
        CFLAGS += -Wextra
endif

all: genetIC

%.o: %.cpp ; $(CXX) $(CFLAGS) $(CODEOPTIONS) $(GIT_VARIABLES) -I$(CPATH) $(FFTW) -c $< -o $@

genetIC: src/main.o src/tools/filesystem.o src/tools/progress/progress.o src/tools/logging.o
		$(CXX) $(CFLAGS) -o genetIC $(GIT_VARIABLES) -I$(CPATH) $(FFTW) src/main.o src/tools/filesystem.o src/tools/progress/progress.o src/tools/logging.o -L$(LPATH) $(GSLFLAGS) -lm $(FFTWLIB)

clean:
	rm -f genetIC
	rm -f src/*.o
	rm -f src/*/*.o
	rm -f src/*/*/*.o
