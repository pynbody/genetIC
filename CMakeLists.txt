cmake_minimum_required(VERSION 3.2)
project(IC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        IC/src/argsort.hpp
        IC/src/constraintapplicator.hpp
        IC/src/cosmo.hpp
        IC/src/dummyic.hpp
        IC/src/fft.hpp
        IC/src/float_types.hpp
        IC/src/grid.hpp
        IC/src/HDF_IO.hh
        IC/src/ic.hpp
        IC/src/io.hpp
        IC/src/main.cpp
        IC/src/mapper.hpp
        IC/src/numpy.hpp
        IC/src/parser.hpp
        IC/src/vectormath.hpp
        IC/src/signaling.hpp
        IC/src/randomfieldgenerator.hpp
        IC/src/filesystem.h
        IC/src/filter.hpp
        IC/src/multilevelcontext.hpp
        IC/src/onelevelconstraint.hpp
        IC/src/filesystem.cpp
        IC/src/multilevelfield.hpp
        IC/src/particle.hpp
        IC/src/progress/progress.cpp
        IC/Makefile)

add_definitions(-DFFTW3 -DFFTW_THREADS)

add_executable(IC ${SOURCE_FILES} )