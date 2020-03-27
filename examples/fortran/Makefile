################################################################################
# Makefile to build the sample executables for the SCARF3D library
################################################################################

include ../../Makefile.inc

# list Fortran objects
OBJS = driver.o scarflib_aux.o

# default target: link to create "driver.exe"
fortran: $(OBJS)
	$(FC) $(FFLAGS) -o driver.exe $(OBJS) $(LIB_LINK_FLAGS)

# list dependencies
driver.o : scarflib_aux.o

# compile objects
$(OBJS) : %.o : %.f90
	$(FC) $(FFLAGS) $(PP_FLAGS) $(INCL) -c -o $@ $<

clean:
	-rm -f *.mod *.o driver.exe fft_* spec_*

# flags
INCL = -I../../include

LIB_LINK_FLAGS  = -L../../lib -lscarf3d
LIB_LINK_FLAGS += -L$(TRNG_DIR)/lib -ltrng4 -lstdc++
LIB_LINK_FLAGS += -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f

# prepare preprocessor directives
PP_FLAGS =

ifeq ($(PRECISION),double)
	PP_FLAGS += -DDOUBLE_PREC
endif

ifeq ($(SPECTRAL),yes)
	PP_FLAGS += -DSPECTRAL
endif

# optimizationflags for each compiler
ifeq ($(COMPILER),gcc)
  FFLAGS   = -O3 -funroll-loops -fopenacc -foffload="-O3 -lm" -foffload=nvptx-none -fopt-info-optimized-omp -cpp
else ifeq ($(COMPILER),pgi)
  FFLAGS   = -fast -acc -ta=tesla:cc75 -cpp
else ifeq ($(COMPILER),intel)
  FFLAGS   = -O3 -unroll-aggressive -qopt-prefetch -cpp
endif

# debugging flags
ifeq ($(DEBUG),yes)
  ifeq ($(COMPILER),gcc)
    FFLAGS   = -Og -fno-lto -g -fcheck=all -fbacktrace -fimplicit-none -Wall -cpp
  else ifeq ($(COMPILER),pgi)
    FFLAGS   = -O0 -g -traceback -Minfo -cpp
  else ifeq ($(COMPILER),intel)
    FFLAGS   = -O0 -g -traceback -check all -check bounds -cpp
  endif
endif
