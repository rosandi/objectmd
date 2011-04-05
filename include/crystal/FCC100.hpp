//--------------FCC 110 Crystal------------------//

#ifndef _CRYSTAL_FCC_100_HPP_
#define _CRYSTAL_FCC_100_HPP_

#include <crystal/Crystalline.hpp>

/**
 @ingroup atom
 @brief FCC(100) Crystal
 
 Direct derivative of Crystal class. This class creates FCC100 surface crystal.
 The algorithm is taken from Impact code. Default=Platinum
*/

class FCC100: public Crystalline {
public:
	FCC100(int XMLayer, int YMLayer, int ZMLayer, string mat_file):
	Crystalline("100", XMLayer, YMLayer, ZMLayer, mat_file) {}

	AtomContainer* Create();
};

#define CrystalFCC100 FCC100

#endif
