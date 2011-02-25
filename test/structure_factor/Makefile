CXX=$(OMD_MPI)/bin/mpic++

EXE=struct

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(OMD_LIBS)

clean:
	rm -rf init.dat info.out  md.out output pomd \
	$(EXE) mem.log run.log *.bin machine
