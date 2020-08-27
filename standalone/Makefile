################################################################################
# Makefile to build the sample executables for the SCARF3D library
################################################################################

include ../Makefile.inc
include ../common_flags.mk

# list Fortran objects
OBJS = main.o scarflib_aux.o

# default target: link to create "driver.exe"
standalone: $(OBJS)
	$(FC) $(FFLAGS) -o scarf3d.exe $(OBJS) $(LIB_LINK_FLAGS)

# list dependencies
main.o : scarflib_aux.o

# compile objects
$(OBJS) : %.o : %.f90
	$(FC) $(FFLAGS) $(PP_FLAGS) $(INCL) -c -o $@ $<

clean:
	-rm -f *.mod *.o scarf3d.exe

# folder with *.mod and *.h files
INCL = -I../include

# path to libraries
LIB_LINK_FLAGS  = -L../lib -l$(LIB_NAME)
LIB_LINK_FLAGS += -L$(TRNG_DIR)/lib -ltrng4 -lstdc++
LIB_LINK_FLAGS += -L$(FFTW_DIR)/lib -lfftw3 -lfftw3f
