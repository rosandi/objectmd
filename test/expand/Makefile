CXX=$(OMD_CXX)
LIB_PATH=-L ../../lib
INCLUDE_PATH=-I ../../include
LIBS=-lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)

EXE=expand

all:
	$(CXX) -O3 -o $(EXE) $(EXE).cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

runsingle:
	./$(EXE)

runmpi: all
	./runmpi

clean:
	rm -rf $(EXE) err me-* system-? true_cry hot cool nohup.* \
	machine  mach-* md.out  pomd Data.??? init_data info.out \
	*.log *.bin
