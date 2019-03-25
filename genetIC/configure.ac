#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([genetIC], [1.0], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE
AC_OUTPUT(Makefile src/Makefile)
AC_CONFIG_SRCDIR([src/ic.hpp])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Options to specify the locations of the required libraries:
# Specify fftw3 directory:
USE_FFTW2="no"
AC_ARG_WITH([fftw3],AS_HELP_STRING([--with-fftw3],[Specify --with-fftw3=/path/to/fftw3 to manually supply the directory location of FFTW3]),[FFTWLIB_DIR=$withval],[])
# Specify fftw2 directory:
AC_ARG_WITH([fftw2],AS_HELP_STRING([--with-fftw2],[Enables the use of the legacy FFTW2 library]),[FFTWLIB_DIR=$withval],[])
# Specify gsl and gslcblas directory:
AC_ARG_WITH([gsl],AS_HELP_STRING([--with-gsl],[Specifies the location of the gsl library as --with-gsl=/path/to/gsl]),[GSL_DIR=$withval USE_FFTW2="yes"],[])

# Checks for libraries.
# Check for the GSL:
AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],[AC_MSG_ERROR([GSLCBLAS not found. Install or specify library location with "--with-gsl=/path/to/lib"])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([GSL not found. Install or specify library location with "--with-gsl=/path/to/lib"])])
# Check for FFTW3, and if not found, check for FFTW2. If neither found, exit with error.
#AS_IF([test $USE_FFTW2="yes"],[AC_CHECK_LIB([fftw3],[fftw_execute],[FFTWLIB="-lfftw3 -lfftw3_threads"],[AC_MSG_ERROR([FFTW3 not found.])])],[AC_CHECK_LIB([dfftw],[fftw_execute],[FFTWLIB="-ldfftw.so -ldrfftw.so"],[AC_MSG_ERROR([FFTW2 not found.])])])
AC_CHECK_LIB([fftw3],[fftw_execute],[FFTWLIB="-lfftw3 -lfftw3_threads"],[AC_CHECK_LIB([dfftw],[fftw_execute],[FFTWLIB="-ldfftw.so -ldrfftw.so"],[AC_MSG_ERROR([No version of FFTW found. Install, or specify directory with --with-fftw3 or --with-fftw2])])])

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h fenv.h float.h inttypes.h limits.h stddef.h stdint.h stdlib.h string.h sys/ioctl.h sys/param.h sys/time.h unistd.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_CHECK_FUNCS([clock_gettime fesetround floor ftruncate getcwd getpagesize gettimeofday localtime_r memmove memset mkdir modf munmap pow select sqrt strcasecmp strerror strstr])


# Substitutes @FFTWLIB@ appearing in Makefile.am with the above specified libraries
AC_SUBST(FFTWLIB)
AC_SUBST(FFTWLIB_DIR)
AC_SUBST(GSL_DIR)

AC_OUTPUT
