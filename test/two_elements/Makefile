CXX=$(OMD_MPI)/bin/mpic++
LIB_PATH=-L ../../lib -L $(OMD_MPI)/lib -L $(OMD_GSL)/lib
INCLUDE_PATH=-I ../../include -I $(OMD_MPI)/include -I $(OMD_GSL)/include
LIBS=-lomd -lomd_object -lm

EXE=two

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

runsingle:
	./$(EXE)

runmpi: all
	./runmpi

clean:
	rm -rf Data.* $(EXE) $(EXE)-p mem.log run.log *.bin info.out md.out \
	machine pomd ar pt init.data
