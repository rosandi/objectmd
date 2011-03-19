CXX=$(OMD_CXX)
LIB_PATH=-L../../lib
INCLUDE_PATH=-I../../include
LIBS=-lm -lomd_object-$(shell uname -m) -lomd-$(shell uname -m)


all: import export prepare pbc

import:
	$(CXX) -g -o import import.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

export:
	$(CXX) -g -o export export.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	
prepare:
	$(CXX) -g -o prepare_crystal prepare_crystal.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	$(CXX) -g -o simtest simtest.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

saveload:
	$(CXX) -g -o save save.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	$(CXX) -g -o load load.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

pbc:
	$(CXX) -g -o pbc_export pbc_export.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)
	$(CXX) -g -o pbc_import pbc_import.cpp $(LIB_PATH) $(INCLUDE_PATH) $(LIBS)

clean:
	rm -rf import export prepare_crystal simtest err me-* system-? true_cry hot cool nohup.* \
	machine  mach-* md.out  pomd Data.??? init_data info.out \
	import.dat imported type0 type1  ID0  ID1 system.cry output.cry mem.log \
	run.log pbc_info.out pbc_init.dat info.out init.dat output \
	hot.load cool.load load save.bin simulation_system.bin pbc_import_restart save pbc_import pbc_export_55.00.bin pbc_export.dat \
	pbc_export pbc_import_p *.dSYM *.info

