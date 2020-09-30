################################################################################
# Makefile to build the sample executables for the SCARF3D library
################################################################################

include ../../Makefile.inc
include ../../common.mk

# list Fortran objects
OBJS = driver.o scarflib_aux.o

# default target: link to create "driver.exe"
fortran: $(OBJS)
	$(FC) $(FFLAGS) -o driverF.exe $(OBJS) $(LIB_LINK_FLAGS)

# list dependencies
driver.o : scarflib_aux.o

# compile objects
$(OBJS) : %.o : %.f90
	$(FC) $(FFLAGS) $(PP_FLAGS) $(INCL) -c -o $@ $<

clean:
	-rm -f *.mod *.o driver.exe fft_* spec_*

# folder with *.mod and *.h files
INCL = -I../../include

# path to libraries
LIB_LINK_FLAGS  = -L../../lib -l$(LIB_NAME)
LIB_LINK_FLAGS += -L$(TRNG_DIR)/lib -ltrng4 -lstdc++
LIB_LINK_FLAGS += -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f
