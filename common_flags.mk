# Users can modify lines below to update compiler flags

# optimization flags for each compiler
ifeq ($(COMPILER),gcc)
  FFLAGS = -O3 -funroll-loops -fopenacc -foffload="-O3 -lm" -foffload=nvptx-none -fopt-info-optimized-omp -cpp
  #CXXFLAGS = -O3 -funroll-loops -cpp
	CFLAGS = -O3 -funroll-loops -fopenacc -foffload="-O3 -lm" -foffload=nvptx-none -fopt-info-optimized-omp -cpp
	CXXFLAGS = -O3 -funroll-loops -fopenacc -foffload="-O3 -lm" -foffload=nvptx-none -fopt-info-optimized-omp -cpp
	LDFLAGS = -lgfortran -lm
else ifeq ($(COMPILER),pgi)
  FFLAGS = -fast -acc -ta=tesla:cc75 -cpp
  #CXXFLAGS = -fast -cpp
	CFLAGS = -fast -acc -ta=tesla:cc75 -cpp
	CXXFLAGS = -fast -acc -ta=tesla:cc75 -cpp
	LDFLAGS = -lpgf90 -lpgf90rtl
else ifeq ($(COMPILER),intel)
  FFLAGS = -O3 -unroll-aggressive -qopt-prefetch -cpp
	CFLAGS = -O3 -unroll-aggressive -qopt-prefetch
  CXXFLAGS = -O3 -unroll-aggressive -qopt-prefetch
	LDFLAGS = -lifcore
endif

# debugging flags
ifeq ($(DEBUG),yes)
  ifeq ($(COMPILER),gcc)
    FFLAGS   = -Og -fno-lto -g -fcheck=all -fbacktrace -fimplicit-none -Wall -cpp
		CFLAGS   = -Og -fno-lto -g -fstack-check -cpp
    CXXFLAGS = -Og -fno-lto -g -fstack-check -cpp
  else ifeq ($(COMPILER),pgi)
    FFLAGS    = -O0 -g -traceback -Minfo -cpp
		CFLAGS    = -O0 -g -Minfo -cpp
    CXXFLAGS  = -O0 -g -Minfo -cpp
  else ifeq ($(COMPILER),intel)
    FFLAGS    = -O0 -g -traceback -check all -check bounds -cpp
		CFLAGS    = -O0 -g -traceback
    CXXFLAGS  = -O0 -g -traceback
  endif
endif

# please do NOT modify below this line!
#-------------------------------------------------------------------------------

# name of library (for linking)
ifeq ($(PRECISION),double)
	LIB_NAME = scarf3d
else
	LIB_NAME = scarf3df
endif

# prepare preprocessor directives
PP_FLAGS =

ifeq ($(PRECISION),double)
	PP_FLAGS += -DDOUBLE_PREC
endif

ifeq ($(SPECTRAL),yes)
	PP_FLAGS += -DSPECTRAL
endif

ifeq ($(DEBUG),yes)
	PP_FLAGS += -DDEBUG
endif

ifeq ($(TIMING),yes)
	PP_FLAGS += -DTIMING
endif

ifeq ($(PFS),yes)
  PP_FLAGS += -DPFS
endif
