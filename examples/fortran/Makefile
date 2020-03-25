################################################################################
# Makefile to build the sample executables for the SCARF3D library
################################################################################

include ../../Makefile.inc

# list Fortran objects
OBJS = driver.o scarflib_aux.o

# default target: link to create "driver.exe"
fortran: $(OBJS)
	$(FC) $(FFLAGS) -o driver.exe $(OBJS) $(LIBS)

# list dependencies
driver.o : scarflib_aux.o

# compile objects
$(OBJS) : %.o : %.f90
	$(FC) $(FFLAGS) $(INCL) -cpp $(CPP) -c -o $@ $<

clean:
	-rm -f *.mod *.o driver.exe fft_* spec_*

# flags
INCL = -I../../include

LIBS = -L../../lib -lscarf3d
LIBS += -L$(TRNG_DIR)/lib -ltrng4 -lstdc++
LIBS += -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f

# prepare preprocessor directives
CPP =

ifeq ($(PRECISION),double)
	CPP += -DDOUBLE_PREC
endif

ifeq ($(SPECTRAL),yes)
	CPP += -DSPECTRAL
endif

# optimizationflags for each compiler
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

# debugging flags
ifeq ($(DEBUG),yes)
  ifeq ($(COMPILER),gcc)
    FFLAGS = -O0 -fno-lto -g -fcheck=all -fbacktrace
    CXXFLAGS = -O0 -fno-lto -g -fbacktrace
  else ifeq ($(COMPILER),pgi)
    FFLAGS    = -O0 -g -traceback -Minfo
    CXXFLAGS  = -g -Minfo
  else ifeq ($(COMPILER),intel)
    FFLAGS    = -O0 -g -traceback -check all -check bounds -debug all
    CXXFLAGS  = -g
  endif
endif
