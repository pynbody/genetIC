cmake_minimum_required(VERSION 3.2)
project(IC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -fopenmp -g -O0")

set(SOURCE_FILES
        genetIC/src/utils.hpp
        genetIC/src/constraints/constraintapplicator.hpp
        genetIC/src/cosmology/parameters.hpp
        genetIC/src/dummyic.hpp
        genetIC/src/numerics/fourier.hpp
        genetIC/src/float_types.hpp
        genetIC/src/grid.hpp
        genetIC/src/ic.hpp
        genetIC/src/io.hpp
        genetIC/src/io/grafic.hpp
        genetIC/src/io/gadget.hpp
        genetIC/src/io/input.hpp
        genetIC/src/io/tipsy.hpp
        genetIC/src/main.cpp
        genetIC/src/coordinate.hpp
        genetIC/src/cosmology/camb.hpp
        genetIC/src/particles/mapper.hpp
        genetIC/src/io/numpy.hpp
        genetIC/src/parser.hpp
        genetIC/src/numerics/vectormath.hpp
        genetIC/src/util/signaling.hpp
        genetIC/src/field/randomfieldgenerator.hpp
        genetIC/src/util/filesystem.h
        genetIC/src/filter.hpp
        genetIC/src/multilevelcontext.hpp
        genetIC/src/constraints/onelevelconstraintgenerator.hpp
        genetIC/src/util/filesystem.cpp
        genetIC/src/field/multilevelfield.hpp
        genetIC/src/particles/particle.hpp
        genetIC/src/progress/progress.cpp
        genetIC/src/field/field.hpp
        genetIC/src/constraints/multilevelconstraintgenerator.hpp
        genetIC/src/numerics/interpolation.hpp
        genetIC/Makefile genetIC/src/particles/zeldovich.hpp genetIC/src/particles/generator.hpp genetIC/src/particles/multilevelgenerator.hpp)

include_directories( /opt/local/include )
link_directories(/opt/local/lib )
include_directories( genetIC/  )
link_libraries(fftw3 m fftw3 fftw3_threads gsl gslcblas)

add_definitions(-DFFTW3 -DFFTW_THREADS -DOPENMP -DGIT_VERSION -DGIT_MODIFIED -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION)
add_compile_options(-Wno-deprecated-declarations)
add_executable(IC ${SOURCE_FILES} )

add_custom_command(TARGET IC
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IC> /Users/app/Science/genetgenetIC/genetIC/IC)