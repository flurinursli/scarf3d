################################################################################
# Configuration file to build the SCARF3D library
################################################################################

# choose COMPILER wrapper
FC = mpif90

# choose COMPILER vendor
COMPILER = pgi

# compile library in single or double precision
PRECISION = single

# compile in debug mode?
DEBUG = no

# activate stopwatch for performance analysis
TIMING = no

# target for parallel filesystem (pfs)?
PFS = no

# enable test of spectral method in driver programs
SPECTRAL = yes

# path to external library TRNG4: update also LD_LIBRARY_PATH if necessary
TRNG_DIR = $(TRNG_PATH)

# path to external library FFTW: update also LD_LIBRARY_PATH if necessary
FFTW_DIR = $(FFTW_PATH)
