################################################################################
# Makefile to build the SCARF3D library
################################################################################

# set dependencies if invoked with 'examples' argument
examples: fortran #c++ c

# invoke MAKE with 'fortran' argument
fortran:
	cd fortran; $(MAKE)	$@

c++:
	cd c++; $(MAKE) $@

c:
	cd c; $(MAKE) $@

clean:
	cd fortran; $(MAKE) $@
	cd c++; $(MAKE) $@
	cd c; $(MAKE) $@

.PHONY: fortran c++ c
