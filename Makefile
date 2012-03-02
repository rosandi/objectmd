# Project: objectmd
# (c) 2005, Yudi Rosandi

include config.make

M_OMD_CXX=$(shell ./select-compiler)

all: omd gadget
#app

debug:
	rm -f lib/*
	$(MAKE) CXX="$(M_OMD_CXX) -g $(OMDFLAG)" -C src/omd
	$(MAKE) CXX="$(M_OMD_CXX) -g $(OMDFLAG)" -C src/gadget

optimize:
	rm -f lib/*
	$(MAKE) CXX="$(M_OMD_CXX) -O3 -Wall -Wno-strict-aliasing $(OMDFLAG)" -C src/omd
	$(MAKE) CXX="$(M_OMD_CXX) -O3 -Wall -Wno-strict-aliasing $(OMDFLAG)" -C src/gadget

omd: src/omd/*.cpp
	rm -f lib/libomd.a
	$(MAKE) CXX="$(M_OMD_CXX) $(OMDFLAG)" -C src/omd

gadget: src/gadget/*.cpp
	rm -f lib/libomd_object.a
	$(MAKE) CXX="$(M_OMD_CXX) $(OMDFLAG)" -C src/gadget
	
app: src/app/*.cpp
	$(MAKE) CXX="$(M_OMD_CXX) $(OMDFLAG)" -C src/app

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
