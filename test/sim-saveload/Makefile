CXX=$(OMD_CXX)
LIB_PATH=-L../../lib
INCLUDE_PATH=-I../../include
LIBS=-lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)


all:
	$(CXX) -g -o pbc_export pbc_export.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	$(CXX) -g -o pbc_import pbc_import.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

clean:
	rm -rf *.dSYM pbc_export pbc_import *.bin *.out *.info

