CXX=$(OMD_CXX)
LIB_PATH=-L../../lib -L/usr/local/lib
INCLUDE_PATH=-I../../include -I/usr/local/include
LIBS=-lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)

EXE=nrb

all:
	$(CXX) -O3 -o $(EXE) $(EXE).cpp $(FLAGS) $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

runsingle:
	./$(EXE)

clean:
	rm -rf Data.* init.* $(EXE) pomd md.out output info.out *.log mach par.md

