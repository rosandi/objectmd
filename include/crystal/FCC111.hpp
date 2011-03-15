//----------------------Crystal FCC(111)----------------------//

#ifndef _CRYSTAL_FCC_111_HPP_
#define _CRYSTAL_FCC_111_HPP_

#include <omd/container.hpp>

/**
  @ingroup atom
  @brief FCC(111) Crystal
  
  The class to create a bulk FCC(111) crystal.  

  Geometry:
 
   <pre>
   xy:                        xz: 
             north   x1,y1               top      x1,z1
          +---------+                +-----------+
          |         |                |           |
   west   |   (.)   | east      west |    (x)    | east
          |   +z    |                |    +y     |
          +---------+                +-----------+
      x0,y0  south               x0,z0  bottom
   </pre>
*/

class CrystalFCC111:public AtomContainer {
protected:
	OMD_INT xml, yml, zml;
	OMD_FLOAT lattice_constant;
	OMD_FLOAT XMLDist, YMLDist, ZMLDist;

public: 
	CrystalFCC111(OMD_INT XMLayer, OMD_INT YMLayer, OMD_INT ZMLayer, string mat_file);
	CrystalFCC111(const OMD_CHAR* fname, string mat_file);
	virtual AtomContainer* Create();
	virtual SysBox& CalcBox();
};


#endif
