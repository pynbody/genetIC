cmake_minimum_required(VERSION 3.2)
project(IC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -fopenmp -g -O0")

set(SOURCE_FILES
        IC/src/utils.hpp
        IC/src/constraints/constraintapplicator.hpp
        IC/src/cosmology/parameters.hpp
        IC/src/dummyic.hpp
        IC/src/numerics/fourier.hpp
        IC/src/float_types.hpp
        IC/src/grid.hpp
        IC/src/ic.hpp
        IC/src/io.hpp
        IC/src/io/grafic.hpp
        IC/src/io/gadget.hpp
        IC/src/io/input.hpp
        IC/src/io/tipsy.hpp
        IC/src/main.cpp
        IC/src/coordinate.hpp
        IC/src/cosmology/camb.hpp
        IC/src/particles/mapper.hpp
        IC/src/io/numpy.hpp
        IC/src/parser.hpp
        IC/src/numerics/vectormath.hpp
        IC/src/util/signaling.hpp
        IC/src/field/randomfieldgenerator.hpp
        IC/src/util/filesystem.h
        IC/src/filter.hpp
        IC/src/multilevelcontext.hpp
        IC/src/constraints/onelevelconstraintgenerator.hpp
        IC/src/util/filesystem.cpp
        IC/src/field/multilevelfield.hpp
        IC/src/particles/particle.hpp
        IC/src/progress/progress.cpp
        IC/src/field/field.hpp
        IC/src/constraints/multilevelconstraintgenerator.hpp
        IC/src/numerics/interpolation.hpp
        IC/Makefile IC/src/particles/zeldovich.hpp IC/src/particles/generator.hpp IC/src/particles/multilevelgenerator.hpp)

include_directories( /opt/local/include )
link_directories(/opt/local/lib )
include_directories( IC/  )
link_libraries(fftw3 m fftw3 fftw3_threads gsl gslcblas)

add_definitions(-DFFTW3 -DFFTW_THREADS -DOPENMP -DGIT_VERSION -DGIT_MODIFIED -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION)
add_compile_options(-Wno-deprecated-declarations)
add_executable(IC ${SOURCE_FILES} )

add_custom_command(TARGET IC
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IC> /Users/app/Science/genetIC/IC/IC)