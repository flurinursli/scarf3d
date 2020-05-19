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
PRECISION = double

# compile in debug mode?
DEBUG = yes

# activate stopwatch for performance analysis
TIMING = yes

# target for parallel filesystem (pfs)?
PFS = no

# enable test of spectral method in driver programs
SPECTRAL = no

# path to external library TRNG4: update also LD_LIBRARY_PATH if necessary
TRNG_DIR = $(TRNG_PATH)

# path to external library FFTW: update also LD_LIBRARY_PATH if necessary
FFTW_DIR = $(FFTW_PATH)
