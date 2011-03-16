CXX=$(OMD_CXX)
EXE=example
LIBS=-L../../lib -I../../include -lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)

all:
	$(CXX) -O3 -o $(EXE) $(EXE).cpp $(LIBS)

clean:
	rm -rf init.dat info.out  md.out output pomd \
	$(EXE) mem.log run.log *.bin machine *.dSYM
