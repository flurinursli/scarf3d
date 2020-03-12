################################################################################
# Makefile to build the SCARF3D library
################################################################################

SCARF3D_DIR = $(CURDIR)

.PHONY : lib examples clean install_dir

# "make" or "make all" to build library and sample programs
all : lib examples

# "make lib" to build library only
lib:
	cd src; $(MAKE) $@

# "make examples"
examples: lib
	cd $@; $(MAKE) $@

clean:
	cd src; $(MAKE) $@
	cd lib; $(MAKE) $@
	cd include; rm -f *.mod
	cd examples; $(MAKE) $@

install_dir:
	mkdir -p $(DESTDIR)$(prefix)
	mkdir -p $(DESTDIR)$(prefix)/include
	mkdri -p $(DESTDIR)$(prefix)/lib

install: all install_dir
	cp $(SCARF3D_DIR)/include/*.mod $(DESTDIR)$(prefix)/include
	cp $(SCARF3D_DIR)/lib/lib*.a $(DESTDIR)$(prefix)/lib
