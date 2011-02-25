CXX=$(OMD_MPI)/bin/mpic++
LIB_PATH=-L ../../lib -L $(OMD_MPI)/lib -L $(OMD_GSL)/lib
INCLUDE_PATH=-I ../../include -I $(OMD_MPI)/include -I $(OMD_GSL)/include
LIBS=-lomd -lomd_object -lm

EXE=sigact

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

runsingle:
	./$(EXE)

runmpi: all
	./runmpi

clean:
	rm -rf data.init  info.out  md.out output  pomd \
	$(EXE) mem.log run.log *.bin machine
