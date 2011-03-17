# Project: objectmd
# (c) 2005, Yudi Rosandi

include config.make

all: omd gadget
#app

debug:
	rm -f lib/*
	$(MAKE) CXX="$(M_OMD_CXX) -g" -C src/omd
	$(MAKE) CXX="$(M_OMD_CXX) -g" -C src/gadget

optimize:
	rm -f lib/*
	$(MAKE) CXX="$(M_OMD_CXX) -O3 -Wall -Wno-strict-aliasing" -C src/omd
	$(MAKE) CXX="$(M_OMD_CXX) -O3 -Wall -Wno-strict-aliasing" -C src/gadget

omd: src/omd/*.cpp
	rm -f lib/libomd.a
	$(MAKE) CXX=$(M_OMD_CXX) -C src/omd

gadget: src/gadget/*.cpp
	rm -f lib/libomd_object.a
	$(MAKE) CXX=$(M_OMD_CXX) -C src/gadget
	
app: src/app/*.cpp
	$(MAKE) CXX=$(M_OMD_CXX) -C src/app

install: optimize
	./install.sh \
	--with-mpi=$(MPI_HOME) \
	--prefix=$(INSTALL_DIR) \
	--conf=$(CONF_FILE)

docu:
	doxygen doc/dox.cfg

doc: docu 

clean:
	$(MAKE) clean -C src/omd
	rm -rf doc/html
