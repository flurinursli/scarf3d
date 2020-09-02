################################################################################
# Configuration file to build the SCARF3D library
################################################################################

# MPI compiler wrappers
FC  = mpif90
C++ = mpic++
CC  = mpicc

# determine compiler suite.
# set manually if variable not defined (default accepted values are "gcc", "pgi", "intel")
COMPILER = $(VENDOR)

# compile library in single (recommended) or double precision,
PRECISION = single

# compile in debug mode? (yes/no)
DEBUG = no

# activate stopwatch for performance analysis? (yes/no)
TIMING = yes

# target for I/O parallel file-system (yes/no)?
PFS = yes

# testing of SRM in driver programs? (yes/no)
SPECTRAL = yes

# allow arbitrarily oriented heterogeneity in FIM? (yes/no)
DIP = no

# path to external library TRNG4.
# set manually if ENV variable not defined (update LD_LIBRARY_PATH if necessary)
TRNG_DIR = $(TRNG_PATH)

# path to external library FFTW.
# set manually if ENV variable not defined (update LD_LIBRARY_PATH if necessary)
FFTW_DIR = $(FFTW_PATH)
