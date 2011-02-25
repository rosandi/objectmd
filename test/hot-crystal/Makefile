CXX=$(OMD_MPI)/bin/mpic++
LIB_PATH=-L ../../lib -L $(HOME)/mpi/lib -L $(HOME)/gsl/lib
INCLUDE_PATH=-I ../../include -I $(HOME)/mpi/include -I $(HOME)/gsl/include
LIBS=-lomd -lomd_object -lm

EXE=hot

all:
#	svn update
	$(CXX) -g -o $(EXE) $(EXE).cpp $(FLAGS) $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
#	strip $(EXE)

runmpi: all
	./runmpi

clean:
	rm -rf $(EXE) err ar pt nohup.* \
	machine  md.out  pomd Data.??? init_data info.out 	 
