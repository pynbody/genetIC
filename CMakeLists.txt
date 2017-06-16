cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_COMPILER /opt/local/bin/g++-mp6)

project(IC)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/genetIC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -fopenmp -g ")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(SOURCE_FILES
        genetIC/src/tools/util_functions.hpp
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
        genetIC/src/simulation/particles/mapper/mapper.hpp
        genetIC/src/io/numpy.hpp
        genetIC/src/tools/parser.hpp
        genetIC/src/tools/numerics/vectormath.hpp
        genetIC/src/tools/signaling.hpp
        genetIC/src/simulation/field/randomfieldgenerator.hpp
        genetIC/src/tools/filesystem.h
        genetIC/src/simulation/filters/filter.hpp
        genetIC/src/simulation/multilevelcontext/multilevelcontext.hpp
        genetIC/src/tools/filesystem.cpp
        genetIC/src/simulation/field/multilevelfield.hpp
        genetIC/src/simulation/particles/particle.hpp
        genetIC/src/tools/progress/progress.cpp
        genetIC/src/simulation/field/field.hpp
        genetIC/src/tools/numerics/interpolation.hpp
        genetIC/Makefile
        genetIC/src/simulation/particles/zeldovich.hpp
        genetIC/src/simulation/particles/generator.hpp
        genetIC/src/simulation/particles/multilevelgenerator.hpp
        genetIC/src/tools/data_types/complex.hpp
        genetIC/src/simulation/filters/filterfamily.hpp
        genetIC/src/simulation/grid/virtualgrid.hpp
        genetIC/src/simulation/particles/mapper/mapperiterator.hpp
        genetIC/src/simulation/particles/mapper/onelevelmapper.hpp
        genetIC/src/simulation/particles/mapper/twolevelmapper.hpp
        genetIC/src/simulation/particles/mapper/gasmapper.hpp
        genetIC/src/simulation/particles/mapper/graficmapper.hpp
        genetIC/src/simulation/field/evaluator.hpp
        genetIC/src/simulation/window.hpp
        genetIC/src/simulation/modifications/modification.hpp
        genetIC/src/simulation/modifications/modificationmanager.hpp
        genetIC/src/simulation/modifications/linearmodification.hpp
        genetIC/src/simulation/modifications/quadraticmodification.hpp)

include_directories( /opt/local/include )
link_directories(/opt/local/lib )
include_directories( genetIC/  )
link_libraries(fftw3 m fftw3 fftw3_threads gsl gslcblas)

add_definitions(-DFFTW3 -DFFTW_THREADS -DOPENMP -DGIT_VERSION -DGIT_MODIFIED -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION)
add_compile_options(-Wextra)
add_executable(IC ${SOURCE_FILES})

add_executable(IC_CFFT ${SOURCE_FILES} )

target_compile_definitions(IC_CFFT PUBLIC -DUSE_COMPLEX_FFT)