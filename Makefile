################################################################################
# Makefile to build the SCARF3D library
################################################################################

SCARF3D_DIR = $(CURDIR)

# set default installation directory (overridden if defined by user)
prefix ?= ~/scarf3d_test

# "make" or "make all" to build library and sample programs
all : lib examples

# "make exe" to build library and standalone application
exe: lib standalone

# "make lib" to build library only
lib:
	cd src; $(MAKE) $@

# "make examples" to build examples (and library, if not yet done)
examples: lib
	cd $@; $(MAKE) $@

# "make standalone" to build standalone application (and library, if not yet done)
standalone: lib
	cd $@; $(MAKE) $@

clean:
	cd src; $(MAKE) $@
	cd examples; $(MAKE) $@
	cd standalone; $(MAKE) $@

# "make install" copy the library to the desired location
install: lib
	mkdir -p $(prefix)
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	cp $(SCARF3D_DIR)/include/*.mod $(prefix)/include
	cp $(SCARF3D_DIR)/include/*.h $(prefix)/include
	cp $(SCARF3D_DIR)/lib/lib*.a $(prefix)/lib

.PHONY : lib examples exe standalone clean install
