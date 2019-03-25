#Adding the git header information:
GIT_VERSION = $(shell git rev-parse --short HEAD)
GIT_MODIFIED = $(shell git ls-files -m | tr "\n" " ")
GIT_VARIABLES = -DGIT_VERSION='"$(GIT_VERSION)"' -DGIT_MODIFIED='"$(GIT_MODIFIED)"'
#Library locations
AM_LDFLAGS = -L@FFTWLIB_DIR@ -L@GSL_DIR@
#Binary to compile:
bin_PROGRAMS = IC
#Source files to compile:
IC_SOURCES = main.cpp tools/progress/progress.cpp tools/filesystem.cpp
IC_CXXFLAGS= -Wall -g -lpthread -Og -fopenmp -std=c++14   -fdiagnostics-color=auto -I`pwd` -DOPENMP  -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -I$(top_srcdir)/ $(GIT_VARIABLES)
IC_LDADD = -lgsl -lgslcblas @FFTWLIB@
