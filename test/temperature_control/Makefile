CXX=$(OMD_MPI)/bin/mpic++
LIB_PATH=-L ../../lib -L $(HOME)/mpi/lib -L $(HOME)/gsl/lib
INCLUDE_PATH=-I ../../include -I $(HOME)/mpi/include -I $(HOME)/gsl/include
LIBS=-lomd -lomd_object -lm

EXE=tempctr

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(FLAGS) $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

runmpi: all
	./runmpi

clean:
	rm -rf $(EXE) machine  md.out  pomd Data.??? init.dat info.out *.log
