//----------------------Crystal FCC----------------------//

#ifndef _CRYSTAL_FCC_HPP_
#define _CRYSTAL_FCC_HPP_

#include <crystal/Crystalline.hpp>

class FCC:public Crystalline {

public: 
	FCC(string orientation, int XMLayer, int YMLayer, int ZMLayer, string mat_file):
	Crystalline(orientation, XMLayer, YMLayer, ZMLayer, mat_file){
		set_name("FCC");
		register_class(get_name());
	}	
	
	void Create111();
	void Create100();
	
	virtual AtomContainer* Create() {

		if(orientation=="111") {
			Create111();
			created=true;
			CalcBox();
			return this;
		}
		
		if(orientation=="100") {
			Create100();
			created=true;
			CalcBox();
			return this;
		}
		
		die("orientation:"+orientation+" is not implemented...");

	}
};

#endif
