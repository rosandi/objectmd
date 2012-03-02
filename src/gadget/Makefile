# Project: objectmd
# (c) 2005, Yudi Rosandi

#CXX = g++ -m64 -O3

include ../../config.make

CXX = g++ -g -Wall
CXXFLAGS = -I../../include $(MPI_FLAGS)
AR = ar
RANLIB = ranlib
RM = rm
MV = mv

LIBOBJ = FCC.o \
	team.o \
	StructureDetector.o \
	VerletList.o \
	VerletListFull.o \
	sw.o


LIBOMD = ../../lib/libomd_object-$(shell uname -m).a

all: $(LIBOMD)

$(LIBOMD): $(LIBOBJ)
	$(AR) -r $(LIBOMD) $(LIBOBJ)
	$(RANLIB) $(LIBOMD)
	$(RM) *.o

clean: 
	${RM} -f *.o ../../lib/$(LIBOMD)

.c.o: 
	$(CXX) -c $*.c $(CXXFLAGS)
