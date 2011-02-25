//--------------FCC 110 Crystal------------------//

#ifndef _CRYSTAL_FCC_100_HPP_
#define _CRYSTAL_FCC_100_HPP_

#include <omd/container.hpp>
#include <crystal/FCC111.hpp>

/**
 @ingroup atom
 @brief FCC(100) Crystal
 
 Direct derivative of Crystal class. This class creates FCC100 surface crystal.
 The algorithm is taken from Impact code. Default=Platinum
*/

class CrystalFCC100: public CrystalFCC111 {
public:
	CrystalFCC100(OMD_INT XMLayer, OMD_INT YMLayer, OMD_INT ZMLayer, string mat_file):
	CrystalFCC111(XMLayer, YMLayer, ZMLayer, mat_file) {}

	/**
	 This part of code is taken from impact code, AG Urbassek (Thomas J. Colla).
	 The function takes monolayer numbers, x, y, z, and the lattice constant.
	 The boundary box of the crystal is not defined here.
	 The caller object must take responsible for it. 
	*/

	AtomContainer* Create();
};

#endif
