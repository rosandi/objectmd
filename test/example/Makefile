CXX=mpicxx
EXE=example
LIBS=-L../../lib -I../../include -lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m) -lmuparser

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(LIBS)

clean:
	rm -rf init.dat info.out  md.out output pomd \
	$(EXE) mem.log run.log *.bin machine *.dSYM
