//----------------------Crystal FCC(111)----------------------//

#ifndef _CRYSTAL_FCC_111_HPP_
#define _CRYSTAL_FCC_111_HPP_

#include <crystal/Crystalline.hpp>

class FCC111:public Crystalline {

public: 
	FCC111(int XMLayer, int YMLayer, int ZMLayer, string mat_file):
	Crystalline("111", XMLayer, YMLayer, ZMLayer, mat_file){}	
	
	virtual AtomContainer* Create();
};

// compatibility...
#define CrystalFCC111 FCC111

#endif
