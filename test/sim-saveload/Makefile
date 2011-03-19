CXX=$(OMD_CXX)
LIB_PATH=-L../../lib
INCLUDE_PATH=-I../../include
LIBS=-lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)


all:
	$(CXX) -g -o export export.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	$(CXX) -g -o import import.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

clean:
	rm -rf *.dSYM export import *.bin *.out *.cry *.info \
	pomd *.log
	

