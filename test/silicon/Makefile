CXX=$(OMD_CXX)
LIBS=-L../../lib -I../../include -lm -lomd-$(shell uname -m) -lomd_object-$(shell uname -m)

all: freerun quench

freerun:
	$(CXX) -g -o freerun freerun.cpp $(LIBS)

quench:
	$(CXX) -g -o quench quench.cpp $(LIBS)

damp:
	$(CXX) -g -o damp damp.cpp $(LIBS)

clean:
	rm -rf init.dat info.out  si.out output pomd quench \
	$(EXE) mem.log run.log *.bin machine *.dSYM Data.* dump-*
