################################################################################
# Makefile for building the SCARF3D library
################################################################################

# Edit here below
#-------------------------------------------------------------------------------
# build
# choose COMPILER wrapper
FC = mpif90
# choose COMPILER vendor
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
TRNG_DIR = $(TRNG_PATH)
# path to external library FFTW: update also LD_LIBRARY_PATH if necessary
FFTW_DIR = $(FFTW_PATH)
#-------------------------------------------------------------------------------
# Do NOT edit rest of file

# prepare paths to "include"
INCLUDE = -I$(TRNG_DIR)/include
INCLUDE += -I$(FFTW_DIR)/include

# prepare libraries flag
LIBRARY = -L$(TRNG_DIR)/lib -ltrng4 -lstdc++
LIBRARY += -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f

# prepare preprocessor directives
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
ifeq ($(debug),yes)
  ifeq ($(COMPILER),gcc)
    FFLAGS = -O0 -fno-lto -g -fcheck=all -fbacktrace
    CXXFLAGS = -O0 -fno-lto -g -fbacktrace
  else ifeq ($(comiler),pgi)
    FFLAGS    = -O0 -g -c -traceback -Minfo
    CXXFLAGS  = -g -Minfo
  else ifeq ($(COMPILER),intel)
    FFLAGS    = -O0 -g -traceback -check all -check bounds -debug all
    CXXFLAGS  = -g
  endif
endif

PROGRAM = scarf3d.exe
MAIN_OBJ = driver.o

# set default target
default : $(PROGRAM)

# list modules
MOD_OBJ = scarflib.o scarflib_common.o scarflib_fft.o scarflib_spec.o scarflib_aux.o

# list of all Fortran object files
F90_OBJ = $(MOD_OBJ)
F90_OBJ += $(MAIN_OBJ)

# compile but do not link
$(F90_OBJ) : %.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -cpp $(CPP) -c -o $@ $<

# list of all C object files
CPP_OBJ = prng.o

# compile but do not link
$(CPP_OBJ) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -cpp $(CPP) -c -o $@ $<

OBJECTS = $(F90_OBJ)
OBJECTS += $(CPP_OBJ)

# now link
# $(PROGRAM) : $(OBJECTS)
#  	$(FC) $(FFLAGS) $(LIBRARY) -o $@ $^
$(PROGRAM) : $(OBJECTS)
	$(FC) $(FFLAGS) -o $(PROGRAM) $(OBJECTS) $(LIBRARY)


# make debug
debug:
	@echo "EXECUTABLE = $(EXECUTABLE)"
	@echo "MAIN_OBJ = $(MAIN_OBJ)"
	@echo "MOD_OBJ = $(MOD_OBJ)"
	@echo "OBJECTS = $(OBJECTS)"
	@echo "LIBRARY = $(LIBRARY)"

clean:
	rm  scarflib.mod scarflib_spec.mod scarflib_fft.mod scarflib_common.mod scarflib_aux.mod
	rm  $(OBJECTS)

.PHONY: debug default clean

# list dependencies for main program
#$(MAIN_OBJ) : $(MOD_OBJ)
$(MAIN_OBJ) : scarflib.o scarflib_common.o scarflib_aux.o

# inter-module dependencies
scarflib_spec.o scarflib_fft.o : scarflib_common.o
scarflib.o: scarflib_common.o scarflib_spec.o scarflib_fft.o
scarflib_aux.o : scarflib_common.o

#
#
# OBJECTS = scarflib_common.o scarflib_fft.o scarflib_spec.o driver.o prng.o
#
# # executable
# scarf3d: $(OBJECTS)
# 	$(FC) $(FFLAGS) -o scarf3d $(OBJECTS) -L${subst :, -L,$(LIBRARY_PATH)} -lstdc++ -ltrng4 -lfftw3f -lfftw3
#
# # mod files
# scarflib.mod: scarflib_common.mod scarflib_fft.mod scarflib_spec.mod scarflib.f90
# 	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib.f90
#
# scarflib_common.mod: scarflib_common.f90
# 	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_common.f90
#
# scarflib_fft.mod: scarflib_common.mod scarflib_fft.f90
# 	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_fft.f90 -I${subst :, -I,$(CPATH)}
#
# scarflib_spec.mod: scarflib_common.mod scarflib_spec.f90
# 	$(FC) $(FFLAGS) -c scarflib_spec.f90
#
# scarflib.o: scarflib_common.mod scarflib_spec.mod scarflib_fft.mod scarflib.f90
# 	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib.f90
#
# scarflib_common.o: scarflib_common.f90
# 	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_common.f90
#
# scarflib_fft.o: scarflib_common.mod scarflib_fft.f90
# 	$(FC) $(FFLAGS) -cpp $(CPP) -c scarflib_fft.f90 -I${subst :, -I,$(CPATH)}
#
# scarflib_spec.o: scarflib_common.mod scarflib_spec.f90
# 	$(FC) $(FFLAGS) -c scarflib_spec.f90
#
# driver.o: scarflib_common.mod scarflib_fft.mod scarflib_spec.mod driver.f90
# 	$(FC) $(FFLAGS) -c driver.f90
#
# # driver.o: scarflib.mod driver.f90
# # 	$(FC) $(FFLAGS) -c driver.f90
#
# #%.o: %.f
# #	$(FC) $(FFLAGS) -c $<
#
# prng.o: prng.cpp
# 	$(CXX) $(CXXFLAGS) -cpp $(CPP) -c prng.cpp
#
#
# clean:
# 	rm  scarflib_common.mod scarflib_spec.mod scarflib_fft.mod scarflib.mod
# 	rm  $(OBJECTS)
