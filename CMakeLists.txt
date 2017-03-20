cmake_minimum_required(VERSION 3.2)
project(IC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -fopenmp -g -O0")

set(SOURCE_FILES
        genetIC/src/tools/util_functions.hpp
        genetIC/src/simulation/constraints/constraintapplicator.hpp
        genetIC/src/cosmology/parameters.hpp
        genetIC/src/dummyic.hpp
        genetIC/src/tools/numerics/fourier.hpp
        genetIC/src/tools/data_types/float_types.hpp
        genetIC/src/simulation/grid/grid.hpp
        genetIC/src/ic.hpp
        genetIC/src/io.hpp
        genetIC/src/io/grafic.hpp
        genetIC/src/io/gadget.hpp
        genetIC/src/io/input.hpp
        genetIC/src/io/tipsy.hpp
        genetIC/src/main.cpp
        genetIC/src/simulation/coordinate.hpp
        genetIC/src/cosmology/camb.hpp
        genetIC/src/simulation/particles/mapper.hpp
        genetIC/src/io/numpy.hpp
        genetIC/src/tools/parser.hpp
        genetIC/src/tools/numerics/vectormath.hpp
        genetIC/src/tools/signaling.hpp
        genetIC/src/simulation/field/randomfieldgenerator.hpp
        genetIC/src/tools/filesystem.h
        genetIC/src/simulation/filters/filter.hpp
        genetIC/src/simulation/multilevelcontext/multilevelcontext.hpp
        genetIC/src/simulation/constraints/onelevelconstraintgenerator.hpp
        genetIC/src/tools/filesystem.cpp
        genetIC/src/simulation/field/multilevelfield.hpp
        genetIC/src/simulation/particles/particle.hpp
        genetIC/src/tools/progress/progress.cpp
        genetIC/src/simulation/field/field.hpp
        genetIC/src/simulation/constraints/multilevelconstraintgenerator.hpp
        genetIC/src/tools/numerics/interpolation.hpp
        genetIC/Makefile
        genetIC/src/simulation/particles/zeldovich.hpp
        genetIC/src/simulation/particles/generator.hpp
        genetIC/src/simulation/particles/multilevelgenerator.hpp
        genetIC/src/tools/data_types/complex.hpp genetIC/src/simulation/filters/filterfamily.hpp)

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