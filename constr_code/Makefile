CC    	= g++
C	= gcc
#mpiCC
CFLAGS	= -Wall -O2
#-ffast-math:
#Sets -fno-math-errno, -funsafe-math-optimizations, -ffinite-math-only, -fno-rounding-math, -fno-signaling-nans and -fcx-limited-range.
#This option causes the preprocessor macro __FAST_MATH__ to be defined.
#This option is not turned on by any -O option besides -Ofast since it can result in incorrect output for programs that depend 
#on an exact implementation of IEEE or ISO rules/specifications for math functions. It may, however, yield faster code for programs 
#that do not require the guarantees of these specifications. 
CPATH = /usr/include
CPATH2= /users/nroth/localcode/fftw/include
LFLAGS  =
#LPATH   = /usr/lib
LPATH   =/usr/local/lib
LPATH2   = /users/nroth/localcode/fftw/lib

#FFTW_FLAGS= 
FFTW_FLAGS= -ldrfftw -ldfftw

#all: ICgd ICg ICgbs ICgbd ICconst
#all: la_test mat_test mat_test_s 
#all: bns_test
#all: cov_iter IC_iter ft3 ICgbd
#all: IC_iter ICgbd IC_iter_c 
#all: mat_test_c IC_iter_c IC_iter_c_single IC_iter_c_ms IC_deriv
all: IC_deriv_int


IC_deriv_int: ICgauss_deriv_clean_int.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_deriv_int -I$(CPATH) ICgauss_deriv_clean_int.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -DHAVE_HDF5
ft3: fftw3_test.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ft3 -I$(CPATH) fftw3_test.cc -L$(LPATH) -lm -lfftw3
cov_iter: cov_iter.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ci -I$(CPATH) cov_iter.cc -L$(LPATH) -lm -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
IC_deriv: ICgauss_deriv_spam.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_deriv -I$(CPATH) ICgauss_deriv_spam.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -DHAVE_HDF5
IC_deriv_alpha: ICgauss_iter_deriv_alpha_test.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_deriv_a -I$(CPATH) ICgauss_iter_deriv_alpha_test.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -DHAVE_HDF5
IC_iter: ICgauss_iter.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_i -I$(CPATH) ICgauss_iter.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
IC_iter_c: ICgauss_iter_clean.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_ic -I$(CPATH) ICgauss_iter_clean.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -ldrfftw -ldfftw -lgsl -lgslcblas -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -DHAVE_HDF5
IC_iter_c_ms: ICgauss_iter_clean_man_shift.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_ic_ms -I$(CPATH) ICgauss_iter_clean_man_shift.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION -DHAVE_HDF5
IC_iter_c_single: ICgauss_iter_clean.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC_ic_single -I$(CPATH) ICgauss_iter_clean.cc -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -lsrfftw -lsfftw -lgsl -lgslcblas -lhdf5 -D H5_USE_16_API -DHAVE_HDF5
bns_test: bns_test.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o bt -I$(CPATH) bns_test.cc -L$(LPATH) -lm -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
la_test: la_test.c Makefile
	cc -O2 -o lat -I/usr/include la_test.c -L/usr/local/lib -lm -framework accelerate
mat_test: mat_test.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o mt -I$(CPATH) mat_test.cc -L$(LPATH) -lm -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
mat_test_c: mat_test_complex.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o mtc -I$(CPATH) mat_test_complex.cc -L$(LPATH) -lm -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
mat_test_s: mat_test.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o mts -I$(CPATH) mat_test.cc -L$(LPATH) -lm -lsrfftw -lsfftw -lgsl -lgslcblas -framework accelerate -lhdf5 -D H5_USE_16_API		
ICconst: ICgauss_both_const.cc Makefile
	$(CC) $(CFLAGS) -o ICconst -I$(CPATH) ICgauss_both_const.cc -L$(LPATH) -lm -ldrfftw -ldfftw -lgsl -lgslcblas -framework accelerate -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
ICgbs: ICgauss_both.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICgbs -I$(CPATH) -L$(LPATH) -lm -lgsl -lgslcblas -I$(CPATH) ICgauss_both.cc -L$(LPATH) -lsrfftw -lsfftw -lhdf5 -D H5_USE_16_API
ICgbd: ICgauss_both.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICgbd -I$(CPATH) -L$(LPATH) -lm -lgsl -lgslcblas -I$(CPATH) ICgauss_both.cc -L$(LPATH) -ldrfftw -ldfftw -lhdf5 -D H5_USE_16_API -DDOUBLEPRECISION -DOUTPUT_IN_DOUBLEPRECISION
ICg: ICgauss.c Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICg -I$(CPATH) -L$(LPATH) -lm -lgsl -lgslcblas -I$(CPATH) ICgauss.c -L$(LPATH) -lsrfftw -lsfftw -lhdf5 -D H5_USE_16_API
ICgd: ICgauss_double.c Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICgd -I$(CPATH) -L$(LPATH) -lm -lgsl -lgslcblas -I$(CPATH) ICgauss_double.c -L$(LPATH) -ldrfftw -ldfftw -lhdf5 -D H5_USE_16_API
IC: IC.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o IC -I$(CPATH) IC.cc -L$(LPATH) -lrfftw -lfftw -lm -lhdf5 -D H5_USE_16_API -lgsl -lgslcblas
ICs: ICs.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICs -I$(CPATH) -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -lgsl -lgslcblas -I$(CPATH) ICs.cc -L$(LPATH) -lsrfftw -lsfftw
ICsp: ICsp.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICsp -I$(CPATH) -L$(LPATH) -lm -lgsl -lgslcblas -I$(CPATH) ICsp.cc -L$(LPATH) -lsrfftw -lsfftw
ICd: ICd.cc Makefile HDF_IO.hh 
	$(CC) $(CFLAGS) -o ICd ICd.cc -I$(CPATH2) -L$(LPATH2) -Wl,-rpath,/users/nroth/localcode/fftw/lib -ldfftw -ldrfftw -I$(CPATH) -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -lgsl -lgslcblas
ICsdnt: ICsdnt.cc Makefile HDF_IO.hh
	$(CC) $(CFLAGS) -o ICsdnt -I$(CPATH) -L$(LPATH) -lm -lhdf5 -D H5_USE_16_API -lgsl -lgslcblas -I$(CPATH2) ICsdnt.cc -L$(LPATH2) -lsrfftw -lsfftw
clean: 
	rm -rf ICs ICsp ICsdnt ICgd ICg ICgbs ICgbd IC_ic IC_ic_ms IC_ic_single
