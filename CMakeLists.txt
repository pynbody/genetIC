cmake_minimum_required(VERSION 3.2)

project(IC)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/genetIC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14  -g ")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -DDEBUG_INFO")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp -O3 -DOPENMP")

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
        genetIC/src/simulation/multilevelgrid/multilevelgrid.hpp
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
        genetIC/src/simulation/particles/offsetgenerator.hpp
        genetIC/src/simulation/field/evaluator.hpp
        genetIC/src/simulation/window.hpp
        genetIC/src/simulation/modifications/modification.hpp
        genetIC/src/simulation/modifications/modificationmanager.hpp
        genetIC/src/simulation/modifications/linearmodification.hpp
        genetIC/src/simulation/modifications/quadraticmodification.hpp
        genetIC/src/simulation/multilevelgrid/mask.hpp
        genetIC/src/tools/memmap.hpp
        genetIC/src/tools/numerics/tricubic.hpp genetIC/src/tools/logging.hpp genetIC/src/tools/logging.cpp genetIC/src/simulation/modifications/splice.hpp)

include_directories( /opt/local/include )
link_directories(/opt/local/lib )
include_directories( genetIC/  )
link_libraries(fftw3 m fftw3 fftw3_threads gsl gslcblas)


exec_program(
        "git"
        ${CMAKE_CURRENT_SOURCE_DIR}
        ARGS "rev-parse --short HEAD"
        OUTPUT_VARIABLE GIT_VERSION
)

exec_program(
        "bash"
        ${CMAKE_CURRENT_SOURCE_DIR}
        ARGS "-c 'git ls-files -m | tr \"\\n\" \" \"'"
        OUTPUT_VARIABLE GIT_MODIFIED
)

# Here are possible code options. Probably don't fiddle with these unless you know what you're doing.
#
# -DCUBIC_INTERPOLATION: Use cubic interpolation when moving information from a coarse to a fine grid. Slower but
#    much better.
# -DVELOCITY_MODIFICATION_GRADIENT_FOURIER_SPACE: As the name suggests, use "ik" as the gradient operator; otherwise
#    uses a fourth order stencil in real space. Enabling this is essential if you disable CUBIC_INTERPOLATION and/or
#    FILTER_ON_COARSE_GRID.
# -DZELDOVICH_GRADIENT_FOURIER_SPACE: Same as VELOCITY_MODIFICATION_GRADIENT_FOURIER_SPACE, but for computing the
#    Zeldovich displacement field.
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
# -DDEBUG_INFO: Output extra debug/performance information
add_definitions(-DCUBIC_INTERPOLATION
        -DFRACTIONAL_K_SPLIT=0.3
        -DFRACTIONAL_FILTER_TEMPERATURE=0.1
        -DFILTER_ON_COARSE_GRID
        -DFFTW3
        -DFFTW_THREADS
        -DGIT_VERSION="${GIT_VERSION}"
        -DGIT_MODIFIED="${GIT_MODIFIED}"
        -DDOUBLEPRECISION
        -DOUTPUT_IN_DOUBLEPRECISION
        -DZELDOVICH_GRADIENT_FOURIER_SPACE)
add_compile_options(-Wextra)
add_executable(genetIC ${SOURCE_FILES})