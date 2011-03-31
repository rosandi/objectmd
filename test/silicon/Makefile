CXX=$(OMD_CXX)
EXE=freerun
LIBS=-L../../lib -I../../include -lm -lomd-$(shell uname -m) -lomd_object-$(shell uname -m)

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(LIBS)

run:
	./run
	
clean:
	rm -rf init.dat info.out  si.out output pomd \
	$(EXE) mem.log run.log *.bin machine *.dSYM Data.* dump-*
