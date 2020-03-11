################################################################################
# Makefile for building the SCARF3D library
################################################################################

# Edit here below
#-------------------------------------------------------------------------------
# build
# choose compiler wrapper
FC = mpif90
# choose compiler vendor
COMPILER = pgi
# compile library in single or double precision
PRECISION = single
# compile in debug mode?
DEBUG = no
# activate stopwatch for performance analysis
TIMING = yes
# target for parallel filesystem (pfs)?
PFS = no
# path to external library TRNG4: update also LD_LIBRARY_PATH if necessary
TRNG_PATH += $(TRNG_PATH)
# path to external library FFTW: update also LD_LIBRARY_PATH if necessary
FFTW_PATH += $(FFTW_PATH)
#-------------------------------------------------------------------------------
# Do NOT edit rest of file

CPP =

ifeq ($(PRECISION),double)
   CPP += -DDOUBLE_PREC
endif

ifeq ($(TIMING),yes)
   CPP += -DTIMING
endif

ifeq ($(PFS),yes)
   CPP += -DPFS
endif

ifeq ($(COMPILER),gcc)
   FFLAGS = -O3 -fopenacc -foffload=-lm
   CXXFLAGS = -O3
else ifeq ($(COMPILER),pgi)
   FFLAGS = -fast -acc -ta=tesla:cc75
   CXXFLAGS = -fast
else ifeq ($(COMPILER),intel)
   FFLAGS = -O3
   CXXFLAGS = -O3
endif

ifeq ($(DEBUG),yes)
   ifeq ($(COMPILER),gcc)
      FFLAGS = -O0 -fno-lto -g -fcheck=all -fbacktrace
      CXXFLAGS = -O0 -fno-lto -g -fbacktrace
   else ifeq ($(COMPILER),pgi)
      FFLAGS    = -O0 -g -c -traceback -Minfo
      CXXFLAGS  = -g -Minfo
   else ifeq ($(COMPILER),intel)
      FFLAGS    = -O0 -g -traceback -check all -check bounds -debug all
      CXXFLAGS  = -g
   endif
endif

OBJECTS = scarflib.o scarflib_common.o prng.o scarflib_fft.o scarflib_spec.o driver.o

# executable
scarf3d: $(OBJECTS)
	$(FC) $(FFLAGS) -o scarf3d $(OBJECTS) -L${subst :, -L,$(LIBRARY_PATH)} -lstdc++ -ltrng4 -lfftw3f -lfftw3

# mod files
scarflib.mod: scarflib_common.mod scarflib_fft.mod scarflib_spec.mod scarflib.f90
	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib.f90

scarflib_common.mod: scarflib_common.f90
	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_common.f90

scarflib_fft.mod: scarflib_common.mod scarflib_fft.f90
	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_fft.f90 -I${subst :, -I,$(CPATH)}

scarflib_spec.mod: scarflib_common.mod scarflib_spec.f90
	$(FC) $(FFLAGS) -c scarflib_spec.f90

scarflib.o: scarflib_common.mod scarflib_spec.mod scarflib_fft.mod scarflib.f90
	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib.f90

scarflib_common.o: scarflib_common.f90
	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_common.f90

scarflib_fft.o: scarflib_common.mod scarflib_fft.f90
	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_fft.f90 -I${subst :, -I,$(CPATH)}

scarflib_spec.o: scarflib_common.mod scarflib_spec.f90
	$(FC) $(FFLAGS) -c scarflib_spec.f90

#driver.o: scarflib_common.mod scarflib_fft.mod scarflib_spec.mod driver.f90
#	$(FC) $(FFLAGS) -c driver.f90

driver.o: scarflib.mod driver.f90
	$(FC) $(FFLAGS) -c driver.f90

#%.o: %.f
#	$(FC) $(FFLAGS) -c $<

prng.o: prng.cpp
	$(CXX) $(CXXFLAGS) -cpp $(CPP) -c prng.cpp


clean:
	rm  scarflib_common.mod scarflib_spec.mod scarflib_fft.mod scarflib.mod
	rm  $(OBJECTS)
