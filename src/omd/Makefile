# Project: objectmd
# (c) 2005, Yudi Rosandi

#CXX = g++ -m64 -O3

include ../../config.make

CXX = g++ -g -Wall
CXXFLAGS = -I../../include -I$(MPI_HOME)/include
AR = ar
RANLIB = ranlib
RM = rm
MV = mv

LIBOBJ = toolkit.o \
	tablereader.o \
	base.o \
	container.o \
	atomgroup.o \
	gadget.o \
	integrator.o \
	comhandler.o \
	system.o \
	systemgrid.o

LIBOMD = ../../lib/libomd-$(shell uname -m).a


all: $(LIBOMD)

$(LIBOMD): $(LIBOBJ)
	$(AR) -r $(LIBOMD) $(LIBOBJ)
	$(RANLIB) $(LIBOMD)
	$(RM) *.o

clean: 
	${RM} -f *.o ../../lib/*

.c.o: 
	$(CXX) -c $*.c $(CXXFLAGS)
