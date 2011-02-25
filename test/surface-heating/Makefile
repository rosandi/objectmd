CXX=$(OMD_MPI)/bin/mpic++
LIB_PATH=-L ../../lib -L $(HOME)/mpi/lib
INCLUDE_PATH=-I ../../include -I $(HOME)/mpi/include
LIBS=-lm -lomd -lomd_object

all:
#	svn update
	$(CXX) -o nonref -O3 nonref.cpp $(FLAGS) $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	strip nonref

clean:
	rm -rf nonref md.out ERR OUT MAKE.OUT output pomd init_data mach nohup.* err
	