CXX=$(OMD_CXX)
LIB_PATH=-L../../lib
INCLUDE_PATH=-I../../include
LIBS=-lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)
EXE=lop

all:
	$(CXX) -g -o $(EXE) $(EXE).cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

test: all
	sh do_test
	        
clean:
	rm -rf init.dat info.out  md.out output pomd \
	$(EXE) mem.log run.log *.bin machine 0K.lop 300K.lop
