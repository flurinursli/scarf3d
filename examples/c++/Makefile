################################################################################
# Makefile to build the sample executables for the SCARF3D library
################################################################################

include ../../Makefile.inc
include ../../common.mk

# list Fortran objects
F90_OBJS = scarflib_aux.o

# compile F90 objects
$(F90_OBJS) : %.o : %.f90
	$(FC) $(FFLAGS) $(INCL) -c -o $@ $<

# list C++ objects
CPP_OBJS = driver.o

# compile C++ objects
$(CPP_OBJS) : %.o : %.cpp
	$(C++) $(CXXFLAGS) $(PP_FLAGS) $(INCL) -c -o $@ $<

# define list of Fortran + C objects
OBJS = $(F90_OBJS)
OBJS += $(CPP_OBJS)

MPI_LINK_FLAGS = $(shell $(FC) --showme:link)

# default target: link to create "driver.exe"
c++: $(OBJS)
#	$(FC) $(FFLAGS) -o driver.exe $(OBJS) $(LIB_LINK_FLAGS) $(LDFLAGS)
	$(C++) $(CXXFLAGS) $(PP_FLAGS) -o driverCPP.exe $(OBJS) $(LIB_LINK_FLAGS) $(MPI_LINK_FLAGS) $(LDFLAGS)


# list dependencies
driver.o : scarflib_aux.o

clean:
	-rm -f *.mod *.o driver.exe fft_* spec_*

# flags
INCL = -I../../include

LIB_LINK_FLAGS  = -L../../lib -l$(LIB_NAME)
LIB_LINK_FLAGS += -L$(TRNG_DIR)/lib -ltrng4 -lstdc++
LIB_LINK_FLAGS += -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f
