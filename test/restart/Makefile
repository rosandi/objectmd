CXX=$(OMD_MPI)/bin/mpic++
INCLUDE_PATH=-I ../../include -I $(OMD_MPI)/include
LIBS=-L../../lib -I../../include -lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)

EXE=small

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp  $(LIBS)
	$(CXX) -g -o $(EXE)-p $(EXE)-p.cpp  $(LIBS)

runsingle:
	./$(EXE)

runmpi: all
	./runmpi

clean:
	rm -rf Data.* $(EXE) $(EXE)-p mem.log run.log save.bin *.bin.* \
	info.out md.out machine pomd nohup.out data.*
