################################################################################
# Makefile to build the sample executables for the SCARF3D library
################################################################################

include ../../Makefile.inc

ifeq ($(PRECISION),double)
	LIB_NAME = scarf3d
else
	LIB_NAME = scarf3df
endif

# list Fortran objects
F90_OBJS = scarflib_aux.o

# compile F90 objects
$(F90_OBJS) : %.o : %.f90
	$(FC) $(FFLAGS) $(INCL) -c -o $@ $<

# list C objects
C_OBJS = driver.o

# compile C objects
$(C_OBJS) : %.o : %.c
	$(CC) $(CFLAGS) $(PP_FLAGS) $(INCL) -c -o $@ $<

# define list of Fortran + C objects
OBJS = $(F90_OBJS)
OBJS += $(C_OBJS)

MPI_LINK_FLAGS = $(shell $(FC) --showme:link)

# default target: link to create "driver.exe"
c: $(OBJS)
#	$(FC) $(FFLAGS) -o driver.exe $(OBJS) $(LIB_LINK_FLAGS) $(LDFLAGS)
	$(CC) $(CFLAGS) $(PP_FLAGS) -o driver.exe $(OBJS) $(LIB_LINK_FLAGS) $(MPI_LINK_FLAGS) $(LDFLAGS)


# list dependencies
driver.o : scarflib_aux.o

clean:
	-rm -f *.mod *.o driver.exe fft_* spec_*

# flags
INCL = -I../../include

LIB_LINK_FLAGS  = -L../../lib -l$(LIB_NAME)
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
  FFLAGS = -O3 -funroll-loops -fopenacc -foffload="-O3 -lm" -foffload=nvptx-none -fopt-info-optimized-omp -cpp
  CFLAGS = -O3 -funroll-loops -fopenacc -foffload="-O3 -lm" -foffload=nvptx-none -fopt-info-optimized-omp -cpp
	LDFLAGS = -lgfortran -lm
else ifeq ($(COMPILER),pgi)
  FFLAGS = -fast -acc -ta=tesla:cc75 -cpp
  CFLAGS = -fast -acc -ta=tesla:cc75 -cpp
	LDFLAGS = -lpgf90 -lpgf90rtl
	#LDFLAGS = -Mnomain
else ifeq ($(COMPILER),intel)
  FFLAGS = -O3 -unroll-aggressive -qopt-prefetch -cpp
  CFLAGS = -O3 -unroll-aggressive -qopt-prefetch
	LDFLAGS = -lifcore
	#LDFLAGS = -nofor-main -cxxlib
endif

# debugging flags
ifeq ($(DEBUG),yes)
  ifeq ($(COMPILER),gcc)
    FFLAGS = -Og -fno-lto -g -fcheck=all -fbacktrace -fimplicit-none -Wall -cpp
    CFLAGS = -Og -fno-lto -g -fstack-check -cpp
  else ifeq ($(COMPILER),pgi)
    FFLAGS    = -O0 -g -traceback -Minfo -acc -ta=tesla:cc75 -cpp
    CFLAGS  = -g -Minfo -acc -ta=tesla:cc75 -cpp
  else ifeq ($(COMPILER),intel)
    FFLAGS    = -O0 -g -traceback -check all -check bounds -cpp
    CFLAGS  = -O0 -g -traceback
  endif
endif
