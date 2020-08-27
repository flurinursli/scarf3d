################################################################################
# Configuration file to build the SCARF3D library
################################################################################

# choose COMPILER wrapper
FC  = mpif90
C++ = mpic++
CC  = mpicc

# choose COMPILER vendor
COMPILER = $(VENDOR)

# compile library in single or double precision
PRECISION = single

# compile in debug mode? (yes/no)
DEBUG = no

# activate stopwatch for performance analysis? (yes/no)
TIMING = yes

# target for I/O parallel file-system (yes/no)?
PFS = no

# enable test of spectral method in driver programs? (yes/no)
SPECTRAL = yes

# path to external library TRNG4: update also LD_LIBRARY_PATH if necessary
TRNG_DIR = $(TRNG_PATH)

# path to external library FFTW: update also LD_LIBRARY_PATH if necessary
FFTW_DIR = $(FFTW_PATH)
